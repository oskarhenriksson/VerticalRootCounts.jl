
export lower_bound_of_maximal_positive_steady_state_count,
    lower_bound_of_maximal_positive_root_count,
    lower_bound_of_maximal_positive_root_count_fixed_a_b_h


struct NongenericDirectionError <: Exception
    msg::String
end
Base.showerror(io::IO, e::NongenericDirectionError) = print(io, e.msg)

struct PositiveRootBoundResult
    bound::Int
    a_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}
    b_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}
    method::Symbol
    TropB::Union{TropicalVariety,Nothing}
    TropL::Union{TropicalLinearSpace,Nothing}
    h::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}    
    stable_intersection::Union{StableIntersectionResult,Nothing}
end
function Base.show(io::IO, ::MIME"text/plain", r::PositiveRootBoundResult)
    header = "Result of positive tropical root bound computation"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, " Lower bound on the maximal number of positive roots: ", r.bound)
    if r.method == :degeneracy
        println(io, " Computation method: degeneracy")
    elseif r.method == :stable_intersection
        println(io, " Computation method: stable intersection of binomial and linear parts")
        println(io, " Choice of parameters a: ", "[", join(r.a_spec, ", "), "]")
        println(io, " Choice of constant terms b: ", "[", join(r.b_spec, ", "), "]")
        println(io, " Choice of perturbation h: ", "[", join(r.h, ", "), "]")
    end
end

@doc raw"""
    lower_bound_of_maximal_positive_root_count_fixed_a_b_h(
    F::AugmentedVerticalSystem,
    a_spec::Union{Vector{Int},Vector{QQFieldElem}},
    b_spec::Union{Vector{Int},Vector{QQFieldElem}}, 
    h::Union{Vector{Int},Vector{QQFieldElem}}; 
    TropB::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)

Compute a lower bound of the maximal positive root count for an augmented vertically parametrized system given by 
the matrices `C`, `M` and `L`, given a fixed choice of constant terms `b_spec`, parameters `a_spec` and shift `h` 
of the tropicalized binomial variety.

# Example
```jldoctest
julia> C = matrix(QQ, [1 -1 -1]);

julia> M = matrix(ZZ, [1 0 2; 0 1 1]);

julia> L = matrix(QQ, [1 1]);

julia> F = AugmentedVerticalSystem(C, M, L);

julia> h = [37,97,18];

julia> a = [839, 562, 13];

julia> b = [71];

julia> bound_result = lower_bound_of_maximal_positive_root_count_fixed_a_b_h(F, a, b, h);

julia> bound_result.bound
3

```
"""
function lower_bound_of_maximal_positive_root_count_fixed_a_b_h(
    F::AugmentedVerticalSystem,
    a_spec::Union{Vector{<:Integer},Vector{QQFieldElem}},
    b_spec::Union{Vector{<:Integer},Vector{QQFieldElem}}, 
    h::Union{Vector{<:Integer},Vector{QQFieldElem}}; 
    TropB::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)

    # Minimal presentation
    C_min, M_min = F.C_min, F.M_min

    # Symbolic coefficient matrix for the augmentation part
    Lb = F.Lb 

    # Dimensions
    n, m, d, r = F.n, F.m, F.d, F.r

    # Check that the input is compatbile with the system
    @req length(b_spec) == d "b_spec must have same length as the number of rows of L"
    @req length(h) == r "h must have same length as the number of monomials"
    @req length(a_spec) == m "a_spec must have same length as the number of columns of M"

    # Set up valuated field for tropicalization
    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:r)

    # Tropicalize the binomial part of the modified system
    if isnothing(TropB)
        binomials = vcat([y[i]-prod(x.^M_min[:,i]) for i=1:r], [z[1]-1])
        TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
        verbose && @info "Tropical binomial variety computed"
    end

    @req check_genericity_of_specialization(Lb, b_spec) "b_spec must be generic"
    Lb_spec = evaluate.(Lb, Ref(b_spec))

    @req check_genericity_of_specialization(C_min, a_spec) "a_spec must be generic"
    C_min_spec = evaluate.(C_min, Ref(a_spec))

    if isnothing(TropL)
        linear_part_matrix = block_diagonal_matrix([Lb_spec, C_min_spec])
        kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
        TropL = tropical_linear_space(kernel_matrix)
        verbose && @info "Tropical linear space computed"
    end

    Σ = perturb_and_intersect_if_transversal(TropL, TropB,
                                perturbation=vcat(zeros(Int, n+1), h), with_multiplicities=false)
    if !Σ.is_transverse 
        throw(NongenericDirectionError("The shift of the tropicalized binomial variety not generic; try a different h vector"))
    end

    # Count how many of the tropical points that are positive
    Ilin = ideal(R, C_min_spec*y) + ideal(R, Lb_spec*vcat(x,z))
    normalized_points = (lcm(denominator.(p)) .* p for p in Σ.points)
    bound = count(
        Oscar.is_initial_positive(Ilin, nu, p) 
        for p in normalized_points
    )
    return PositiveRootBoundResult(
        bound, 
        a_spec, 
        b_spec, 
        :stable_intersection, 
        TropB, 
        TropL, 
        h, 
        Σ
    )
end

@doc raw"""
    lower_bound_of_maximal_positive_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix; 
    num_a_b_attempts::Int=5, num_h_attempts_per_a_b::Int=10, verbose::Bool=false)

Computes a lower bound on the maximal positive root count of the augmented vertically parametrized 
system given by the coefficient matrix `C`, the exponent matrix `M`, and the affine form matrix `L`.

The function randomly samples `num_a_b_attempts` choices of the b and k parameters, and
for each such choice `num_h_attempts_per_a_b` shifts of the tropicalized binomial variety 
in the space of auxiliary variables in the modification.


"""
function lower_bound_of_maximal_positive_root_count(F::AugmentedVerticalSystem;
    num_a_b_attempts::Int=5, 
    num_h_attempts_per_a_b::Int=10, 
    show_progress::Bool=true,
    max_entry_size::Int=1000,
    verbose::Bool=false
)

    # Minimal presentation
    C_min, M_min = F.C_min, F.M_min

    # Matrices for the augmentation part
    L = F.L
    Lb = F.Lb 

    # Dimensions
    n, m, d, r = F.n, F.m, F.d, F.r

    # Check whether there are nondegenerate zeros at all
    if !has_nondegenerate_zero(F)
        return PositiveRootBoundResult(
            0, 
            nothing, 
            nothing, 
            :degeneracy,
            nothing, 
            nothing, 
            nothing,
            nothing
        )
    end

    # Tropicalize the binomial part of the modified system
    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:r)
    binomials = vcat([y[i]-prod(x.^M_min[:,i]) for i=1:r], [z[1]-1])
    TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
    verbose && @info "Tropical binomial variety computed"
   
    # Try different choices of b and h
    # Keep track of the maximal positive root count found and associated b and h values
    # Todo: Make this interruptible!
    best_result = nothing

    progress = ProgressMeter.Progress(num_a_b_attempts; 
        dt=0.4, 
        desc="Trying parameter values...", 
        barlen=30,
        output = stdout,
        enabled = show_progress
    );

    for a_b_attempt=1:num_a_b_attempts

        # Pick a generic a
        a_spec = nothing
        while true
            a_spec = rand(1:max_entry_size, m)
            is_generic = check_genericity_of_specialization(C_min, a_spec)
            if is_generic
                break
            end
        end
        C_min_spec = evaluate.(C_min, Ref(a_spec))

        # Pick a generic b
        b_spec = nothing
        while true
            b_spec = L*rand(1:max_entry_size, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
            if is_generic
                break
            end
        end
        Lb_spec = evaluate.(Lb, Ref(b_spec))
    
        # Tropicalize the linear part of the modified system
        linear_part_matrix = block_diagonal_matrix([Lb_spec, C_min_spec])
        kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
        TropL = tropical_linear_space(kernel_matrix)
        verbose && @info "Tropical linear space computed"
    
        # Compute the stable intersection for different h values
        new_result = nothing
        for _ = 1:num_h_attempts_per_a_b
            generic_perturbation = false
            while !generic_perturbation
                try
                    h = rand(1:max_entry_size, r)
                    new_result = lower_bound_of_maximal_positive_root_count_fixed_a_b_h(
                        F, a_spec, b_spec, h; TropB=TropB, TropL=TropL, verbose=verbose
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

            # Update the current best result if new one is better
            if isnothing(best_result) || new_result.bound > best_result.bound
                best_result = new_result
            end
        end

        # Update the progress bar
        ProgressMeter.update!(progress, a_b_attempt; 
            showvalues = [
                ("Number of b attempts", "$(a_b_attempt) ($(num_a_b_attempts))"), 
                ("Current maximal count", best_result.bound)
            ]
        )
    end
    return best_result
end


@doc raw"""
    lower_bound_of_maximal_positive_steady_state_count(rn::ReactionSystem; kwargs...)

    Computes a lower bound on the maximal number of isolated positive steady states 
    that a mass action network `rn` can have.
"""
lower_bound_of_maximal_positive_steady_state_count(rn::ReactionSystem; kwargs...) = 
    lower_bound_of_maximal_positive_root_count(steady_state_system(rn); kwargs...)
