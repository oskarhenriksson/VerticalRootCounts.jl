
export lower_bound_of_maximal_positive_steady_state_count,
    lower_bound_of_maximal_positive_root_count,
    lower_bound_of_maximal_positive_root_count_fixed_b_k_h


struct NongenericDirectionError <: Exception
    msg::String
end
Base.showerror(io::IO, e::NongenericDirectionError) = print(io, e.msg)

struct PositiveRootBound
    bound::Int
    b_spec
    k_spec
    h
end
function Base.show(io::IO, ::MIME"text/plain", r::PositiveRootBound)
    header = "Positive tropical root bound"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, " Lower bound on the maximal number of positive roots: ", r.bound)
    println(io, " Choice of constant terms b: ", r.b_spec)
    println(io, " Choice of parameters k: ", r.k_spec)
    print(io, " Choice of perturbation h: ", r.h)
end

@doc raw"""
    lower_bound_of_maximal_positive_root_count_fixed_b_k_h(
    C::QQMatrix, M::ZZMatrix, L::QQMatrix,
    b_spec::Union{Vector{Int},Vector{QQFieldElem}}, 
    k_spec::Union{Vector{Int},Vector{QQFieldElem}},
    h::Union{Vector{Int},Vector{QQFieldElem}}; 
    TropB::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)

Compute a lower bound of the maximal positive root count for an augmented vertically parametrized system given by 
the matrices `C`, `M` and `L`, given a fixed choice of constant terms `b_spec`, parameters `k_spec` and shift `h` 
of the tropicalized binomial variety.

# Example
```jldoctest
julia> C = matrix(QQ, [1 -1 -1]);

julia> M = matrix(ZZ, [1 0 2; 0 1 1]);

julia> L = matrix(QQ, [1 1]);

julia> h = [37,97,18];

julia> k = [839, 562, 13];

julia> b = [71];

julia> lower_bound_of_maximal_positive_root_count_fixed_b_k_h(C, M, L, b, k, h)
3
```
"""
function lower_bound_of_maximal_positive_root_count_fixed_b_k_h(
    C::QQMatrix, M::ZZMatrix, L::QQMatrix,
    b_spec::Union{Vector{Int},Vector{QQFieldElem}}, 
    k_spec::Union{Vector{Int},Vector{QQFieldElem}},
    h::Union{Vector{Int},Vector{QQFieldElem}}; 
    TropB::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)

    n = nrows(M) #number of variables
    m = ncols(M) #number of parameters
    s = rank(C) #rank
    d = n-s #corank

    C_tilde, M_tilde = minimal_presentation(C, M)
    r = ncols(M_tilde)

    @req nrows(L) == d "L must have the same number of rows as the corank of C"
    @req length(b_spec) == d "b_spec must have same length as the number of rows of L"
    @req length(h) == r "h must have same length as the number of columns of M_tilde"
    @req length(k_spec) == m "k_spec must have same length as the number of columns of M"

    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:r)

    # Tropicalize the binomial part of the modified system
    if isnothing(TropB)
        binomials = vcat([y[i]-prod(x.^M_tilde[:,i]) for i=1:r], [z[1]-1])
        TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
        verbose && @info "Tropical binomial variety computed"
    end

    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(B, d, 1, b))
    @req check_genericity_of_specialization(Lb, b_spec) "b_spec must be generic"
    Lb_spec = evaluate.(Lb, Ref(b_spec))

    @req check_genericity_of_specialization(C_tilde, k_spec) "k_spec must be generic"
    C_tilde_spec = evaluate.(C_tilde, Ref(k_spec))

    if isnothing(TropL)
        linear_part_matrix = block_diagonal_matrix([Lb_spec, C_tilde_spec])
        kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
        TropL = tropical_linear_space(kernel_matrix)
        verbose && @info "Tropical linear space computed"
    end

    result = perturb_and_intersect_if_transversal(TropL, TropB,
                                perturbation=vcat(zeros(Int, n+1), h), with_multiplicities=false)
    if !result.is_transversal 
        throw(NongenericDirectionError("The shift of the tropicalized binomial variety not generic; try a different h vector"))
    end

    # Count how many of the tropical points that are positive
    Ilin = ideal(R, C_tilde_spec*y) + ideal(R, Lb_spec*vcat(x,z))

    normalized_points = (lcm(denominator.(p)) .* p for p in result.points)

    return count(
        Oscar.is_initial_positive(Ilin, nu, p)
        for p in normalized_points
    )
end


@doc raw"""
    lower_bound_of_maximal_positive_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix; 
    num_b_k_attempts::Int=5, num_h_attempts_per_b_k::Int=10, verbose::Bool=false)

Computes a lower bound on the maximal positive root count of the augmented vertically parametrized 
system given by the coefficient matrix `C`, the exponent matrix `M`, and the affine form matrix `L`.

The function randomly samples `num_b_k_attempts` choices of the b and k parameters, and
for each such choice `num_h_attempts_per_b_k` shifts of the tropicalized binomial variety 
in the space of auxiliary variables in the modification.


"""
function lower_bound_of_maximal_positive_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M)); 
    num_b_k_attempts::Int=5, 
    num_h_attempts_per_b_k::Int=10, 
    show_progress::Bool=true,
    max_entry_size::Int=1000,
    verbose::Bool=false
)

    n = nrows(M) #number of variables
    m = ncols(M) #number of parameters
    s = rank(C) #rank
    d = n-s #corank

    C_tilde, M_tilde = minimal_presentation(C, M)
    r = ncols(M_tilde)

    # Check whether there are nondegenerate zeros at all
    if !has_nondegenerate_zero(C, M, L)
        return PositiveRootBound(0, L*rand(1:max_entry_size, n), rand(1:max_entry_size, m), rand(1:max_entry_size, r))
    end

    @req nrows(L) == d "L must have the same number of rows as the corank of C"

    # Tropicalize the binomial part of the modified system
    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:r)
    binomials = vcat([y[i]-prod(x.^M_tilde[:,i]) for i=1:r], [z[1]-1])
    TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
    verbose && @info "Tropical binomial variety computed"
   
    # Try different choices of b and h
    # Keep track of the maximal positive root count found and associated b and h values
    # Todo: Make this interruptible!
    best_count = 0
    best_b = nothing
    best_k = nothing
    best_h = nothing
    progress = ProgressMeter.Progress(num_b_k_attempts; 
        dt=0.4, 
        desc="Trying parameter values...", 
        barlen=30,
        output = stdout,
        enabled = show_progress
    );

    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(B, d, 1, b))
    for b_k_attempt=1:num_b_k_attempts

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

        # Pick a generic k
        k_spec = nothing
        while true
            k_spec = rand(1:max_entry_size, m)
            is_generic = check_genericity_of_specialization(C_tilde, k_spec)
            if is_generic
                break
            end
        end
        C_tilde_spec = evaluate.(C_tilde, Ref(k_spec))
    
        # Tropicalize the linear part of the modified system
        linear_part_matrix = block_diagonal_matrix([Lb_spec, C_tilde_spec])
        kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
        TropL = tropical_linear_space(kernel_matrix)
        verbose && @info "Tropical linear space computed"
    
        # Compute the stable intersection for different h values
        new_count = nothing 
        h = nothing
        for h_attempt = 1:num_h_attempts_per_b_k
            generic_perturbation = false
            while !generic_perturbation
                try
                    h = rand(1:max_entry_size, r)
                    new_count = lower_bound_of_maximal_positive_root_count_fixed_b_k_h(
                        C, M, L, b_spec, k_spec, h; TropB=TropB, TropL=TropL, verbose=verbose
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
            if new_count > best_count || isnothing(best_b) || isnothing(best_k) || isnothing(best_h)
                best_count = new_count
                best_b = b_spec
                best_k = k_spec
                best_h = h
            end
        end

        # Update the progress bar
        ProgressMeter.update!(progress, b_k_attempt; 
            showvalues = [
                ("Number of b attempts", "$(b_k_attempt) ($(num_b_k_attempts))"), 
                ("Current maximal count", best_count)
            ]
        )
    end
    return PositiveRootBound(best_count, best_b, best_k, best_h)
end

@doc raw"""
    lower_bound_of_maximal_positive_steady_state_count(rn::ReactionSystem; kwargs...)

    Computes a lower bound on the maximal number of isolated positive steady states 
    that a mass action network `rn` can have.
"""
lower_bound_of_maximal_positive_steady_state_count(rn::ReactionSystem; kwargs...) = 
    lower_bound_of_maximal_positive_root_count(steady_state_system(rn)...; kwargs...)
