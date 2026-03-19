export steady_state_degree,
    generic_root_count,
    generic_degree

struct GenericRootCountResult
    count::Int
    a_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}
    b_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}
    method::Symbol
    TropB::Union{TropicalVariety,Nothing}
    TropL::Union{TropicalLinearSpace,Nothing}
    h::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}}
    stable_intersection::Union{StableIntersectionResult,Nothing}
end

function Base.show(io::IO, ::MIME"text/plain", r::GenericRootCountResult)
    header = "Result of generic root count computation"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, " Generic root count: ", r.count)
    println(io, " Choice of constant terms b: ", "[", join(r.b_spec, ", "), "]")
    println(io, " Choice of parameters k: ", "[", join(r.a_spec, ", "), "]")
    if r.method == :degeneracy
        println(io, " Computation method: degeneracy")
    elseif r.method == :cotransversality
        println(io, " Computation method: mixed volume")
    elseif r.method == :stable_intersection
        println(io, " Computed method: stable intersection")
        println(io, " Computation of perturbation h: ", "[", join(r.h, ", "), "]")
    end
end



@doc raw"""
generic_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix; 
    a_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}} = nothing,
    b_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}} = nothing, 
    check_transversality::Bool=true, 
    verbose::Bool=false)

Compute the generic root count of an augmented vertical system `F`.

# Example

```jldoctest
julia> C = matrix(QQ, [1 -1 -1]);

julia> M = matrix(ZZ, [1 0 2; 0 1 1]);

julia> L = matrix(QQ, [1 1]);

julia> F = AugmentedVerticalSystem(C, M, L);

julia> generic_root_count(F)
Result of generic root count computation
========================================
 Generic root count: 3
 Choice of constant terms b: [503]
 Choice of parameters k: [802, 980, 905]
 Computation method: mixed volume

julia> generic_root_count(F, check_transversality=false)
Result of generic root count computation
========================================
 Generic root count: 3
 Choice of constant terms b: QQFieldElem[636]
 Choice of parameters k: [854, 593, 711]
 Computation method: stable intersection
 Choice of perturbation h: Int16[-14763, 25201, -30895, 15903, -31687, -7376]
```
"""
function generic_root_count(F::AugmentedVerticalSystem;
        b_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}} = nothing, 
        a_spec::Union{Nothing,Vector{<:Integer},Vector{QQFieldElem}} = nothing,
        check_transversality::Bool=true, 
        verbose::Bool=false)

    @req is_square(F) "The system needs to be square (number of rows of C and L need to sum to the number of varialbes)"

    # Check whether there are nondegenerate zeros at all
    if !has_nondegenerate_zero(F)
        return GenericRootCountResult(
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

    C, M, L = F.C, F.M, F.L #defining matrices
    C_tilde, M_tilde = F.C_tilde, F.M_tilde #minimal presentation
    Lb = F.Lb #symbolic coefficient matrix for the augmentation of the system

    n = F.n #number of variables
    m = F.m #number of parameters
    s = F.s #number of polynomials in the nonlinear part
    d = F.d #number of augmenting linear forms
    r = F.r #number of monomials in the nonlinear part


    # Pick a generic specialization of the parameters
    if isnothing(a_spec)
        is_generic = false
        while !is_generic
            a_spec = rand(1:1000, m)
            is_generic = check_genericity_of_specialization(C_tilde, a_spec)
        end
    end
    @req check_genericity_of_specialization(C_tilde, a_spec) "Choice of parameters needs to be generic"
    C_tilde_spec = evaluate.(C_tilde, Ref(a_spec))

    # Pick a generic specialization of the constant terms
    if isnothing(b_spec)
        is_generic = false
        while !is_generic
            b_spec = L*rand(1:1000, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
        end
    end
    @req check_genericity_of_specialization(Lb, b_spec) "Choice of constant terms needs to be generic"
    Lb_spec = evaluate.(Lb, Ref(b_spec))

    # Compute the generic root count as a mixed volume if the linear part gives a transversal matroid
    if check_transversality
        tp_nonlinear = transversal_presentation(C_tilde_spec)
        tp_affine = transversal_presentation(Lb_spec)
        if tp_nonlinear != false && tp_affine != false
            verbose && @info "Transversal presentations found"
            nonlinear_supports = Matrix{Int}[Matrix{Int}(M_tilde[:,indices]) for indices in tp_nonlinear]
            affine_supports =  Matrix{Int}[hcat([i in 1:n ? standard_vector(i, n) : zeros(Int, n) for i in indices]...) for indices in tp_affine]
            supports = vcat(nonlinear_supports, affine_supports)
            return GenericRootCountResult(
                mixed_volume(supports), 
                a_spec, 
                b_spec, 
                :cotransversality, 
                nothing, 
                nothing,
                nothing, 
                nothing
            )
        end
    end

    # Tropicalize the linear part of the modified system
    linear_part_matrix = block_diagonal_matrix([Lb_spec, C_tilde_spec])
    kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
    TropL = tropical_linear_space(kernel_matrix)
    verbose && @info "Tropical linear space computed"

    # Tropicalize the binomial part of the modified system
    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:r)
    binomials = vcat([y[i]-prod(x.^M_tilde[:,i]) for i=1:r], [z[1]-1])
    TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
    verbose && @info "Tropical binomial variety computed"
 
    # Run perturb_and_intersect_if_transversal until done
    Σ = nothing
    while true
        Σ = perturb_and_intersect_if_transversal(TropL, TropB)
        Σ.is_transverse && break
    end
    return GenericRootCountResult(
        sum(Σ.multiplicities),
        a_spec, 
        b_spec, 
        :stable_intersection,
        TropB,
        TropL,
        Σ.perturbation,
        Σ
    )
end


@doc raw"""
    steady_state_degree(rn::ReactionSystem; kwargs...)

Compute the steady state degree (in the complex torus) of the steady state system of a mass action network `rn`.

# Example
```jldoctest
julia> rn = @reaction_network begin
    k1, X1 --> X2
    k2, X2 --> X1
    k3, 2*X1 + X2 --> 3*X1
end;

julia> steady_state_degree(rn)
Result of generic root count computation
========================================
 Generic root count: 3
 Choice of constant terms b: QQFieldElem[360]
 Choice of parameters k: [530, 377, 969]
 Computation method: mixed volume
````

"""
steady_state_degree(rn::ReactionSystem; kwargs...) = 
    generic_root_count(steady_state_system(rn); kwargs...)




"""
    generic_degree(C::QQMatrix, M::ZZMatrix)

Compute the generic degree of the ideal of a purely vertical system `F`.

In accordance with Bézout's theorem, we compute this by intersecting with a 
generic affine space of complementary dimension.

"""
function generic_degree(F) 

    @req is_purely_vertical(F) "The system needs to be purely vertical (number of rows of L needs to be zero)"

    n = nrows(F.M) #number of variables
    s = rank(F.C) #rank

    # Check for nondegeneracy
    if !has_nondegenerate_zero(F)
        return GenericRootCountResult(0, nothing, nothing, false, nothing, nothing, nothing)
    end

    # Augment the system to a square system by an L with full support
    L_generic = matrix(QQ, rand(Int16, n-s, n))

    # Check that the matroid is uniform (all Plücker coordinates are nonzero)
    all(!is_zero, minors(L_generic, nrows(L_generic)))

    F_augmented = AugmentedVerticalSystem(F.C, F.M, L_generic)

    # Compute the generic root count
    return generic_root_count(F_augmented)

end