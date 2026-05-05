
export toric_root_bound,
    toric_lower_bound_of_maximal_positive_root_count_fixed_b_h,
    toric_lower_bound_of_maximal_positive_root_count

@doc raw"""
ToricRootBoundResult

Result type of the toric root bound computation.

Fields:
- `bound::Int`: The computed upper bound on the number of roots in the complex torus.
- `b_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}`: The choice of constant terms b that achieved the bound (if method is not :degeneracy) 
- `method::Symbol`: The method used to compute the bound. One of :degeneracy, :cotransversality, or :stable_intersection.
- `Trop_toric::Union{TropicalVariety,Nothing}`: The tropical variety of the toric part of the system (if method is :stable_intersection)
- `TropL::Union{TropicalLinearSpace,Nothing}`: The tropical linear space of the linear part of the system (if method is :stable_intersection)
- `h::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}`: The choice of perturbation h that achieved the bound (if method is :stable_intersection)
- `stable_intersection::Union{StableIntersectionResult,Nothing}`: The result of the stable intersection computation (if method is :stable_intersection)
- `cotranversal_presentation_Lb::Union{Nothing,Vector{Vector{Int}}}`: The row supports of the cotransversal presentation for the linear part (if method is :cotransversality)
"""

struct ToricRootBoundResult
    bound::Int
    b_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}
    method::Symbol
    Trop_toric::Union{TropicalVariety,Nothing}
    TropL::Union{TropicalLinearSpace,Nothing}
    h::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}
    stable_intersection::Union{StableIntersectionResult,Nothing}
    cotranversal_presentation_Lb::Union{Nothing,Vector{Vector{Int}}}
end

function Base.show(io::IO, ::MIME"text/plain", r::ToricRootBoundResult)
    header = "Result of toric root bound computation"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, " Toric root bound: ", r.bound)
    if r.method == :degeneracy
        println(io, " Computation method: degeneracy")
        return
    end
    if r.method == :cotransversality 
        println(io, " Computation method: mixed volume for cotransversal presentation")
    elseif r.method == :stable_intersection
            println(io, " Computation method: stable intersection of binomial and linear parts")
    end
    println(io, " Choice of constant terms b: ", "[", join(r.b_spec, ", "), "]")
    if r.method == :stable_intersection
        println(io, " Choice of perturbation h: ", "[", join(r.h, ", "), "]")
    end
    if r.method == :cotransversality
        println(io, " Row supports of cotransversal presentation for the linear part: ")
        for indices in r.cotranversal_presentation_Lb
            println(io, "  [", join(indices, ", "), "]")
        end
    end
end



@doc raw"""
    toric_root_bound(A::ZZMatrix, F::AugmentedVerticalSystem;
    b_spec::Union{Nothing,Vector{Int},Vector{QQFieldElem}}=nothing,
    check_cotransversality::Bool=true,
    verbose::Bool=false)

Given an augmented vertical system `F` that is parametrically toric with 
respect to an exponent matrix A, compute an upper bound on the number of 
roots in the complex torus.

Input:
- `A::ZZMatrix`: Exponent matrix of the toric part of the system (size d×n)
- `F::AugmentedVerticalSystem`: An augmented vertical system with n variables and d augmenting linear forms.

Optional inputs:
- `b_spec::Union{Nothing,Vector{Int},Vector{QQFieldElem}}=nothing`: Optional choice of constant terms for the augmenting linear forms. If not provided, a random generic choice will be made.
- `check_cotransversality::Bool=true`: Whether to check for cotransversality.
- `verbose::Bool=false`: Whether to print detailed information about the computation process.

"""
function toric_root_bound(A::ZZMatrix, F::AugmentedVerticalSystem;
    b_spec::Union{Nothing,Vector{Int},Vector{QQFieldElem}}=nothing,
    check_cotransversality::Bool=true,
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
        return ToricRootBoundResult(0, b_spec, :degeneracy, nothing, nothing, nothing, nothing)
    end

    # Check for transversality
    if check_cotransversality
        tp = cotransversal_presentation(Lb_spec)
        if tp != false
            verbose && @info "Cotransversal presentation found"
            supports = [Matrix{Int}(A_extended[:,indices]) for indices in tp]
            supports_shifted = [S .- min.(0, vec(minimum(S, dims=2))) for S in supports];
            degA = Int(prod(diagonal(snf(A))))
            bound = Int(mixed_volume(supports_shifted)/degA)
            return ToricRootBoundResult(
                bound,
                b_spec,
                :cotransversality,
                nothing,
                nothing,
                nothing,
                nothing,
                tp
            )
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
    Σ = nothing
    while true
        Σ = perturb_and_intersect_if_transversal(TropL, Trop_toric)
        Σ.is_transverse && break
    end

    return ToricRootBoundResult(
        sum(Σ.multiplicities), 
        b_spec, 
        :stable_intersection, 
        Trop_toric, 
        TropL, 
        Σ.perturbation,
        Σ,
        nothing
    )
end


@doc raw"""
PositiveToricRootBoundResult

Result type of the toric lower bound of maximal positive root count. 

Fields:
- `bound::Int`: The computed lower bound on the maximal number of positive roots.
- `b_spec::Union{Vector{QQFieldElem}, Nothing}`: The choice of constant terms b that achieved the bound (if method is :stable_intersection)
- `h::Union{Vector{QQFieldElem}, Nothing}`: The choice of perturbation h that achieved the bound (if method is :stable_intersection)
- `method::Symbol`: The method used to compute the bound.
- `TropB::Union{TropicalVariety, Nothing}`: The tropical variety.
- `TropL::Union{TropicalLinearSpace, Nothing}`: The tropical linear space.
- `stable_intersection::Union{StableIntersectionResult,Nothing}`: The result of the stable intersection computation.

"""
struct PositiveToricRootBoundResult
    bound::Int
    b_spec::Union{Vector{QQFieldElem}, Nothing}
    h::Union{Vector{QQFieldElem}, Nothing}
    method::Symbol
    TropB::Union{TropicalVariety, Nothing}
    TropL::Union{TropicalLinearSpace, Nothing}
    stable_intersection::Union{StableIntersectionResult,Nothing}
end

function Base.show(io::IO, ::MIME"text/plain", r::PositiveToricRootBoundResult)
    header = "Result of positive toric root bound computation"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, " Lower bound on the maximal number of positive roots: ", r.bound)
    if r.method == :degeneracy
        println(io, " Computation method: degeneracy")
    elseif r.method == :stable_intersection
        println(io, " Computation method: stable intersection of binomial and linear parts")
        println(io, " Choice of constant terms b: ", "[", join(r.b_spec, ", "), "]")
        println(io, " Choice of perturbation h: ", "[", join(r.h, ", "), "]")
    end
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

    Σ = perturb_and_intersect_if_transversal(TropL, Trop_toric, perturbation=h, with_multiplicities=false)
    
    if !Σ.is_transverse 
        throw(NongenericDirectionError("Input perturbation not generic"))
    end

    # Count how many of the tropical points that are positive
    Lb_spec = hcat(L, -matrix(QQ.(b_spec)))
    Ilin = ideal(R, Lb_spec*vcat(x,z))
    normalized_points = (lcm(denominator.(p)) .* p for p in Σ.points)
    bound =  count(
        Oscar.is_initial_positive(Ilin, tropical_semiring_map(QQ), p) 
        for p in normalized_points
    )
    return PositiveToricRootBoundResult(bound, b_spec, h, :stable_intersection, Trop_toric, TropL, Σ)
end





@doc raw"""
    toric_lower_bound_of_maximal_positive_root_count(
        A::ZZMatrix, 
        F::AugmentedVerticalSystem; 
        num_b_attempts::Int=5, 
        num_h_attempts_per_b::Int=10,
        target_bound::Union{Nothing,Int}=nothing, 
        max_entry_size::Int=1000,
        show_progress::Bool=true,
        verbose::Bool=false
    )

Given an augmented vertical system `F` that is parametrically toric with 
respect to an exponent matrix A, compute a lower bound on the maximal number 
of positive roots. This is done by trying different choices of the constant 
terms `b` and perturbations `h`. 

Input:
- `A::ZZMatrix`: Exponent matrix of the toric part of the system (size d×n)
- `F::AugmentedVerticalSystem`: An augmented vertical system with n variables and d augmenting linear forms.

Optional inputs:
- `num_b_attempts::Int=5`: Number of different choices of the constant terms `b` to try.
- `num_h_attempts_per_b::Int=10`: Number of different choices of the perturbation `h` to try for each choice of `b`.
- `target_bound::Union{Nothing,Int}=nothing`: If provided, the function will stop as soon as a lower bound greater than or equal to `target_bound` is found.
- `max_entry_size::Int=1000`: Maximum absolute value of the entries of `b` and `h` when they are randomly generated.
- `show_progress::Bool=true`: Whether to show a progress bar during the computation.
- `verbose::Bool=false`: Whether to print detailed information about the computation process.

"""
function toric_lower_bound_of_maximal_positive_root_count(
    A::ZZMatrix, 
    F::AugmentedVerticalSystem; 
    num_b_attempts::Int=5, 
    num_h_attempts_per_b::Int=10,
    target_bound::Union{Nothing,Int}=nothing, 
    max_entry_size::Int=1000,
    show_progress::Bool=true,
    verbose::Bool=false
)

    n, d = F.n, F.d
    L, Lb = F.L, F.Lb
    @req ncols(A) == n "Number of columns of A needs to match the number of variables in the system"
    @req rank(A) == d "System needs to be effectively square"

    if !has_nondegenerate_zero(F)
        return PositiveToricRootBoundResult(
            0, 
            nothing, 
            nothing, 
            :degeneracy,
            nothing, 
            nothing, 
            nothing
        )
    end

    R, x, z = polynomial_ring(QQ, "x"=>1:n, "z"=>1:1)

    # Tropicalize the toric variety
    A_extended = hcat(A, zero_matrix(ZZ,nrows(A),1)) # Add a column corresponding to the homogenization variable
    I_toric = toric_ideal(R, transpose(A_extended))
    Trop_toric = Oscar.tropical_variety_binomial(I_toric,tropical_semiring_map(QQ))
    verbose && @info "Tropicalization of toric variety computed"
  
    # Try different choices of b and h
    # Keep track of the maximal positive root count found and associated b and h values
    # Todo: Make this interruptible!
    best_result = nothing
    progress = ProgressMeter.Progress(num_b_attempts; 
        dt=0.4, 
        desc="Trying parameter values...", 
        barlen=30,
        output = stdout,
        enabled = show_progress
    );

    target_reached = false
    TropL = nothing
    
    for b_attempt=1:num_b_attempts

        if target_reached
            break
        end

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
        if isnothing(TropL)
            TropL = tropical_linear_space(ideal(Lb_spec*vcat(x,z)))
            verbose && @info "Tropical linear space computed"
        end
        
        # Compute the stable intersection for different h values
        new_result = nothing 
        for h_attempt = 1:num_h_attempts_per_b
            generic_perturbation = false
            while !generic_perturbation
                try
                    h = rand(1:max_entry_size, (n+1))
                    new_result = toric_lower_bound_of_maximal_positive_root_count_fixed_b_h(
                        A, F, b_spec, h, Trop_toric=Trop_toric, TropL=TropL, verbose=verbose
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
            if isnothing(best_result) || new_result.bound > best_result.bound
                best_result = new_result
                verbose && @info "New best bound found: $(best_result.bound)\nb = $(best_result.b_spec)\nh = $(best_result.h)"
                if !isnothing(target_bound) && best_result.bound >= target_bound
                    verbose && @info "Target bound reached"
                    target_reached = true
                    break
                end
            end
            
            # Update the progress bar
            ProgressMeter.update!(progress, b_attempt*(num_h_attempts_per_b-1) + h_attempt; 
                showvalues = [
                    ("Number of b attempts", "$(b_attempt) ($(num_b_attempts))"), 
                    ("Current maximal count",  isnothing(best_result) ? 0 : best_result.bound)
                ]
            )
        end

    end
    return best_result
end
