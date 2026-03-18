export steady_state_degree,
    generic_root_count,
    generic_degree
    
@doc raw"""
generic_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix; 
    b_spec = nothing, check_transversality::Bool=true, verbose::Bool=false)

    Compute the generic root count of an augmented vertically parametrized system given by 
    the coefficient matrix `C`, the exponent matrix `M`, and the affine form matrix `L`.

    # Example
    ```jldoctest
    julia> C = matrix(QQ, [1 -1 -1]);

    julia> M = matrix(ZZ, [1 0 2; 0 1 1]);
 
    julia> L = matrix(QQ, [1 1]);

    julia> generic_root_count(C, M, L)
    3

    julia> generic_root_count(C, M, L, check_transversality=false)
    3
    ```

"""
function generic_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M)); 
        b_spec = nothing, 
        a_spec = nothing,
        check_transversality::Bool=true, 
        verbose::Bool=false)

    # Check whether there are nondegenerate zeros at all
    if !has_nondegenerate_zero(C, M, L)
        return 0
    end
    
    n = nrows(M) #number of variables
    m = ncols(M) #number of parameters
    s = rank(C) #rank
    d = nrows(L) #corank

    @req s + d == n "The system needs to be square (number of rows of C and L need to sum to the number of varialbes)"

    # Monomial re-embedding of the system
    C_tilde, M_tilde = minimal_presentation(C, M)
    r = ncols(M_tilde)

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

    # Symbolic coefficient matrix for the augmentation of the system
    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(B, d, 1, b))

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
            return mixed_volume(supports)
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
    rootCountComputed = false
    mults = Int[]
    while !rootCountComputed
        result = perturb_and_intersect_if_transversal(TropL, TropB)
        rootCountComputed = result.is_transversal
        mults = result.multiplicities
    end
    return sum(mults)
end



"""
    generic_root_count(F::AugmentedVerticalSystem; kwargs...)

    Compute the generic root count of an augmented vertical system `F`.

     # Example
    ```jldoctest
    julia> C = matrix(QQ, [1 -1 -1]);

    julia> M = matrix(ZZ, [1 0 2; 0 1 1]);
 
    julia> L = matrix(QQ, [1 1]);

    julia> F = AugmentedVerticalSystem(C, M, L);

    julia> generic_root_count(F)
    3

    julia> generic_root_count(F, check_transversality=false)
    3
    ```

"""
generic_root_count(F::AugmentedVerticalSystem; kwargs...) = generic_root_count(F.C, F.M, F.L; kwargs...)



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
3
````

"""
steady_state_degree(rn::ReactionSystem; kwargs...) = 
    generic_root_count(steady_state_system(rn); kwargs...)




"""
    generic_degree(C::QQMatrix, M::ZZMatrix)

Compute the generic degree of the ideal of a vertical system.

In accordance with Bézout's theorem, we compute this by intersecting with a 
generic affine space of complementary dimension.

"""
function generic_degree(C::QQMatrix, M::ZZMatrix) 

    n = nrows(M) #number of variables
    m = ncols(M) #number of parameters
    s = rank(C) #rank

    # Check for nondegeneracy
    if !has_nondegenerate_zero(C, M)
        return 0
    end

    # Augment the system to a square system by an L with full support
    L_generic = matrix(QQ, rand(Int16, n-s, nrows(M)))

    # Check that the matroid is uniform (all Plücker coordinates are nonzero)
    all(!is_zero, minors(L_generic, nrows(L_generic)))

    # Compute the generic root count
    return generic_root_count(C, M, L_generic)

end

function generic_degree(F::AugmentedVerticalSystem)
    @req nrows(F.L) == 0 "The system needs to be purely vertical (number of rows of L needs to be zero)"
    generic_degree(F.C, F.M)
end
