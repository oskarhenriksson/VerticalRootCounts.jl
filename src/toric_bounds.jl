
export toric_root_bound,
    toric_lower_bound_of_maximal_positive_root_count_fixed_b_h,
    toric_lower_bound_of_maximal_positive_root_count

function toric_root_bound(A::ZZMatrix, L::QQMatrix;
    b_spec::Union{Nothing,Vector{Int},Vector{QQFieldElem}}=nothing,
    check_transversality::Bool=true,
    verbose::Bool=false
)

    n = ncols(A) # number of variables
    d = nrows(L) # number of affine equations
    @req rank(A) == d "System needs to be effectively square"

    # Add a column corresponding to the homogenization variable
    A_extended = hcat(A, zero_matrix(ZZ,nrows(A),1))

    # Pick a generic choice of constant terms
    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(b))

    # Pick a generic specialization of the constant terms
    if isnothing(b_spec)
        is_generic = false
        while !is_generic
            b_spec = L*rand(Int16, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
        end
    end
    @req check_genericity_of_specialization(Lb, b_spec) "Choice of constant terms to be generic"
    Lb_spec = evaluate.(Lb, Ref(b_spec))

    # Nondegeneracy check
    if !has_nondegenerate_zero(Lb_spec, A_extended)
        return 0
    end

    # Check for transversality
    if check_transversality
        tp = transversal_presentation(Lb_spec)
        if tp != false
            verbose && @info "Transversal presentation found"
            supports = [Matrix{Int}(A_extended[:,indices]) for indices in tp]
            supports_shifted = [S .- min.(0, vec(minimum(S, dims=2))) for S in supports];
            degA = Int(prod(diagonal(snf(A))))
            return Int(mixed_volume(supports_shifted)/degA)
        end
    end

    R, x, z = polynomial_ring(QQ, "x"=>1:n, "z"=>1:1)

    # Tropicalize the affine linear space
    TropL = tropical_linear_space(ideal(Lb_spec*vcat(x,z)))
    verbose && @info "Tropical linear space computed"

    # Tropicalize the toric variety
    I_toric = toric_ideal(R, transpose(A_extended))
    Trop_toric = Oscar.tropical_variety_binomial(I_toric,tropical_semiring_map(QQ))
    verbose && @info "Tropicalization of toric variety computed"
   
    # Stable intersection
    rootCountComputed = false
    mults = Int[]
    while !rootCountComputed
        _, rootCountComputed, _, mults = perturb_and_intersect_if_transversal(TropL, Trop_toric)
    end
    return sum(mults)
end


function toric_lower_bound_of_maximal_positive_root_count_fixed_b_h(
    A::ZZMatrix, L::QQMatrix,
    b_spec::Union{Vector{Int},Vector{QQFieldElem}}, 
    h::Union{Vector{Int},Vector{QQFieldElem}}; 
    Trop_toric::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)
    n = ncols(A) # number of variables
    d = nrows(L) # number of affine equations
    @req rank(A) == d "System needs to be effectively square"
    R, x, z = polynomial_ring(QQ, "x"=>1:n, "z"=>1:1)

    if isnothing(TropL)
        # Tropicalize the affine linear space
        Lb_spec = hcat(L, -matrix(QQ.(b_spec)))
        TropL = tropical_linear_space(ideal(Lb_spec*vcat(x,z)))
        verbose && @info "Tropical linear space computed"
    end

    # Tropicalize the toric variety
    if isnothing(Trop_toric)
        A_extended = hcat(A, zero_matrix(ZZ,nrows(A),1)) # Add a column corresponding to the homogenization variable
        I_toric = toric_ideal(R, transpose(A_extended))
        Trop_toric = Oscar.tropical_variety_binomial(I_toric,tropical_semiring_map(QQ))
        verbose && @info "Tropicalization of toric variety computed"
    end

    _, isTransversal, pts, _ = perturb_and_intersect_if_transversal(TropL, Trop_toric, perturbation=h, with_multiplicities=false)
    
    if !isTransversal 
        throw(NongenericDirectionError("Input perturbation not generic"))
    end

    # Count how many of the tropical points that are positive
    Lb_spec = hcat(L, -matrix(QQ.(b_spec)))
    Ilin = ideal(R, Lb_spec*vcat(x,z))
    return count(Oscar.is_initial_positive(Ilin, tropical_semiring_map(QQ), lcm(denominator.(p)) .* p) for p in pts)
end


function toric_lower_bound_of_maximal_positive_root_count(A::ZZMatrix, L::QQMatrix,; 
    num_b_attempts::Int=5, 
    num_h_attempts_per_b::Int=10, 
    show_progress::Bool=true,
    verbose::Bool=false
)
    n = ncols(A) # number of variables
    d = nrows(L) # number of affine equations
    R, x, z = polynomial_ring(QQ, "x"=>1:n, "z"=>1:1)
    @req rank(A) == d "System needs to be effectively square"


    # Tropicalize the toric variety
    A_extended = hcat(A, zero_matrix(ZZ,nrows(A),1)) # Add a column corresponding to the homogenization variable
    I_toric = toric_ideal(R, transpose(A_extended))
    Trop_toric = Oscar.tropical_variety_binomial(I_toric,tropical_semiring_map(QQ))
    verbose && @info "Tropicalization of toric variety computed"
  
    # Try different choices of b and h
    # Keep track of the maximal positive root count found and associated b and h values
    # Todo: Make this interruptible!
    best_count = 0
    best_b = nothing 
    best_h = nothing
    progress = ProgressMeter.Progress(num_b_attempts; 
        dt=0.4, 
        desc="Trying parameter values...", 
        barlen=30,
        output = stdout,
        enabled = show_progress
    );
    for b_attempt=1:num_b_attempts

        # Pick a generic b
        B, b = rational_function_field(QQ, "b"=>1:d)
        Lb = hcat(B.(L), -matrix(B, d, 1, b))
        while true
            b_spec = L*rand(Int16, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
            if is_generic
                Lb_spec = evaluate.(Lb, Ref(b_spec))
                break
            end
        end
        
        # Tropical linear space
        TropL = tropical_linear_space(ideal(Lb_spec*vcat(x,z)))
        verbose && @info "Tropical linear space computed"
    
        # Compute the stable intersection for different h values
        new_count = nothing 
        h = nothing
        for h_attempt = 1:num_h_attempts_per_b
            generic_perturbation = false
            while !generic_perturbation
                try
                    h = rand(1:1000, (n+1))
                    new_count = toric_lower_bound_of_maximal_positive_root_count_fixed_b_h(
                        A, L, b_spec, h, Trop_toric=Trop_toric, TropL=TropL
                    )
                    generic_perturbation = true
                catch err
                    if err isa NongenericDirectionError
                        continue
                    else
                        rethrow()
                    end
                end
            end
              
            # Update the current best count
            if new_count > best_count || isnothing(best_b) || isnothing(best_h)
                best_count = new_count
                best_b = b_spec
                best_h = h
            end
        end

        # Update the progress bar
        ProgressMeter.update!(progress, b_attempt; 
            showvalues = [
                ("Number of b attempts", "$(b_attempt) ($(num_b_attempts))"), 
                ("Current maximal count", best_count)
            ]
        )
    end
    return best_count, best_b, best_h
end

