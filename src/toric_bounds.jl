
export toric_root_bound,
    toric_lower_bound_of_maximal_positive_root_count_fixed_b_h,
    toric_lower_bound_of_maximal_positive_root_count

function toric_root_bound(A::ZZMatrix, F::AugmentedVerticalSystem;
    b_spec::Union{Nothing,Vector{Int},Vector{QQFieldElem}}=nothing,
    check_transversality::Bool=true,
    verbose::Bool=false
)

    L, Lb = F.L, F.Lb
    n, d = F.n, F.d

    @req ncols(A) == n "Number of columns of A needs to match the number of variables in the system"
    @req rank(A) == d "System needs to be effectively square"

    # Add a column corresponding to the homogenization variable
    A_extended = hcat(A, zero_matrix(ZZ,nrows(A),1))

    # Pick a generic specialization of the constant terms
    if isnothing(b_spec)
        is_generic = false
        while !is_generic
            b_spec = L*rand(Int16, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
        end
    end
    @req check_genericity_of_specialization(Lb, b_spec) "Choice of constant terms to be generic"
    @req length(b_spec) == d "b_spec must have same length as the number of rows of L"
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
        result = perturb_and_intersect_if_transversal(TropL, Trop_toric)
        rootCountComputed = result.is_transversal
        mults = result.multiplicities
    end
    return sum(mults)
end



struct ToricPositiveRootBound
    bound::Int
    b_spec::Vector{QQFieldElem}
    h::Vector{QQFieldElem}
    TropB::TropicalVariety
    TropL::TropicalLinearSpace
end

function Base.show(io::IO, ::MIME"text/plain", r::ToricPositiveRootBound)
    header = "Toric positive root bound"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, " Lower bound on the maximal number of positive roots: ", r.bound)
    println(io, " Choice of constant terms b: ", r.b_spec)
    println(io, " Choice of perturbation h: ", r.h)
end


function toric_lower_bound_of_maximal_positive_root_count_fixed_b_h(
    A::ZZMatrix, F::AugmentedVerticalSystem,
    b_spec::Union{Vector{Int},Vector{QQFieldElem}}, 
    h::Union{Vector{Int},Vector{QQFieldElem}}; 
    Trop_toric::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)
    
    n, d = F.n, F.d
    L, Lb = F.L, F.Lb
    @req ncols(A) == n "Number of columns of A needs to match the number of variables in the system"
    @req rank(A) == d "System needs to be effectively square"
    @req length(b_spec) == d "b_spec must have same length as the number of rows of L"
    @req length(h) == n+1 "h must have same length as the number of variables plus one (for the homogenization variable)"   

    R, x, z = polynomial_ring(QQ, "x"=>1:n, "z"=>1:1)

    # Check genericity of input
    @req check_genericity_of_specialization(Lb, b_spec) "b_spec must be generic"

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

    result = perturb_and_intersect_if_transversal(TropL, Trop_toric, perturbation=h, with_multiplicities=false)
    
    if !result.is_transversal 
        throw(NongenericDirectionError("Input perturbation not generic"))
    end

    # Count how many of the tropical points that are positive
    Lb_spec = hcat(L, -matrix(QQ.(b_spec)))
    Ilin = ideal(R, Lb_spec*vcat(x,z))
    normalized_points = (lcm(denominator.(p)) .* p for p in result.points)
    bound =  count(
        Oscar.is_initial_positive(Ilin, tropical_semiring_map(QQ), p) 
        for p in normalized_points
    )
    return ToricPositiveRootBound(bound, b_spec, h, Trop_toric, TropL)
end





function toric_lower_bound_of_maximal_positive_root_count(A::ZZMatrix, F::AugmentedVerticalSystem; 
    num_b_attempts::Int=5, 
    num_h_attempts_per_b::Int=10, 
    max_entry_size::Int=1000,
    show_progress::Bool=true,
    verbose::Bool=false
)

    n, d = F.n, F.d
    L, Lb = F.L, F.Lb
    @req ncols(A) == n "Number of columns of A needs to match the number of variables in the system"
    @req rank(A) == d "System needs to be effectively square"

    R, x, z = polynomial_ring(QQ, "x"=>1:n, "z"=>1:1)

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
    best_TropL = nothing
    progress = ProgressMeter.Progress(num_b_attempts; 
        dt=0.4, 
        desc="Trying parameter values...", 
        barlen=30,
        output = stdout,
        enabled = show_progress
    );
    
    for b_attempt=1:num_b_attempts

        # Pick a generic b
        b_spec = nothing
        while true
            b_spec = L*(rand(1:max_entry_size, n))
            is_generic = check_genericity_of_specialization(Lb, b_spec)
            if is_generic
                break
            end
        end
        Lb_spec = evaluate.(Lb, Ref(b_spec))
        
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
                    h = rand(1:max_entry_size, (n+1))
                    result = toric_lower_bound_of_maximal_positive_root_count_fixed_b_h(
                        A, F, b_spec, h, Trop_toric=Trop_toric, TropL=TropL, verbose=verbose
                    )
                    new_count = result.bound
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
                best_TropL = TropL
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
    return ToricPositiveRootBound(best_count, best_b, best_h, Trop_toric, best_TropL)
end

