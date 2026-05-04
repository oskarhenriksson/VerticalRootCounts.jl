export check_genericity_of_specialization

@doc raw"""
    check_genericity_of_specialization(A::AbstractAlgebra.Generic.MatSpaceElem{<:RingElem}, a_spec::Vector{<:Union{Int, RingElem}})

Check if a symbolic matrix is generic at a given specialization, in the sense that the column matroid of the specialized matrix is the same as the column matroid of the generic matrix.

# Examples
```jldoctest
julia> using Oscar;

julia> K, a = rational_function_field(QQ, "a"=>1:3);

julia> A = matrix(K, [[a[1], a[2], a[3]],[a[1], a[1]+a[2], a[2]+a[3]]]);

julia> check_genericity_of_specialization(A, [1, 1, 1])
false

julia> check_genericity_of_specialization(A, [1, 1, 2])
true
```
"""
function check_genericity_of_specialization(A::AbstractAlgebra.Generic.MatSpaceElem{<:RingElem}, a_spec::Vector{<:Union{Int, RingElem}})::Bool

    R = base_ring(A)
    @req ngens(R) == length(a_spec) "Specialization needs to match the number of parameters"

    # Specialize the symbolic matrix at a_spec
    A_spec = evaluate.(A, Ref(a_spec))

    # Compute its matroid
    M_spec = matroid_from_matrix_columns(A_spec)

    # Check that each circuit of M_spec is dependent in the generic matroid
    for C in circuits(M_spec)
        A_C = A[:,C]
        # Check that all lenght(C) minors of A_C are zero
        if any(!is_zero, minors_iterator(A_C, length(C)))
            return false
        end
    end

    return true
end