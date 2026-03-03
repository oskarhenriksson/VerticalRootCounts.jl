export minimal_presentation, 
    augmented_vertical_system,
    has_nondegenerate_zero

function augmented_vertical_system(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M)))

    KB, k, b = rational_function_field(QQ, "k"=>1:ncols(M), "b"=>1:nrows(L))
    R, x = polynomial_ring(KB, "x"=>1:nrows(M))

    steady_state = C*[k[j]*prod(x .^ M[:, j]) for j in 1:ncols(M)]
    conservation_laws = L*x .- b

    return [steady_state; conservation_laws]
end

"""
    minimal_presentation(C, M)

    Given a vertical system given by`C` and exponent matrix `M`, 
    compute the defining matrices for the monomial re-embedding.

    # Example
    ```jldoctest
    julia> C = matrix(QQ, [0 0 1 -1 1 0; 1 -1 -1 0 0 0; 0 0 0 1 -1 -1])
    [0    0    1   -1    1    0]
    [1   -1   -1    0    0    0]
    [0    0    0    1   -1   -1]

    julia> M = matrix(ZZ, [1 0 0 0 0 0; 0 0 0 1 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 1 1 0 0 0; 0 0 0 0 1 1])
    [1   0   0   0   0   0]
    [0   0   0   1   0   0]
    [1   0   0   0   0   0]
    [0   0   0   1   0   0]
    [0   1   1   0   0   0]
    [0   0   0   0   1   1]

    julia> C_tilde, M_tilde = minimal_presentation(C, M);

    julia> C_tilde
    [   0           k[3]   -k[4]           k[5]]
    [k[1]   -k[2] - k[3]       0              0]
    [   0              0    k[4]   -k[5] - k[6]]

    julia> M_tilde
    [1   0   0   0]
    [0   0   1   0]
    [1   0   0   0]
    [0   0   1   0]
    [0   1   0   0]
    [0   0   0   1]

    ```
"""
function minimal_presentation(C::QQMatrix, M::ZZMatrix)
    columns = [M[:, i] for i in 1:ncols(M)]
    unique_columns = unique(columns)
    r = length(unique_columns)

    M_tilde = matrix(ZZ, hcat(unique_columns...))
    K, k = rational_function_field(QQ, "k"=>1:ncols(M))
    C_tilde = zero_matrix(K, nrows(C), r) 
    for i = 1:r
        indices = findall(c -> c == unique_columns[i], columns)
        for j in indices
            C_tilde[:, i] += k[j] .* C[:, j]
        end
    end
    return C_tilde, M_tilde
end



"""
    has_nondegenerate_zero(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M));
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)

    Check if the augmented vertical system given by `C`, `M`, and `L` (optional input) has a non-degenerate zero.

    This uses the rank condition of [^FHP25].

    # Example

    ```jldoctest
    julia> C = C = matrix(QQ, [[1,-1,-1]]);
    julia> M = matrix(ZZ, [[1,0,2], [0,1,1]]);
    julia> L = matrix(QQ, [[1,1]]);
    julia> has_nondegenerate_zero(C, M, L)
    true
    ```

    ```jldoctest
    julia> C = matrix(QQ, [-1  0  0  0  1  0; 0 -1  0  0  0  1;  0  0 -1  1  0  0]);
    julia> M = matrix(ZZ, [3  2  1  0  0  0; 0  1  0  2  1  0; 0  0  1  0  1  2]);
    julia> has_nondegenerate_zero(C, M)
    false
    ```  
    
    [^FHP25] Feliu, E., Henriksson, O., Pascual-Escudero, B. "Generic consistency and nondegeneracy of vertically parametrized systems". Journal of Algebra 677 (2025).

"""
function has_nondegenerate_zero(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M));
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)
    C = (rref(C)[2])[1:rank(C), :]
    L = rref(L)[2]
    @req ncols(C) == ncols(M) "C and M need to have the same number of columns"
    @req ncols(L) == nrows(M) "L needs to have the same number of columns as M has rows"
    s = rank(C)
    G = kernel(C, side=:right)
    for _ in 1:number_of_attempts
        u = rand(-max_entry_size:max_entry_size, ncols(G))
        nrows(L) == 0 ? h = ones(Int, nrows(M)) : h = rand(-max_entry_size:max_entry_size, nrows(M))
        degeneracy_matrix = vcat(C * diagonal_matrix(G * u) * transpose(M) * diagonal_matrix(h), L)
        if rank(degeneracy_matrix) == s + nrows(L)
            return true
        end
    end
    if certify
        R, u, h = polynomial_ring(QQ, "u" => 1:ncols(G), "h" => 1:nrows(M))
        nrows(L) == 0 ? h = ones(Int, nrows(M)) : nothing
        symbolic_degeneracy_matrix = vcat(C * diagonal_matrix(G * u) * transpose(M) * diagonal_matrix(h), R.(L))
        return !all(is_zero, minors_iterator(symbolic_degeneracy_matrix, s + nrows(L)))
    else
        return false
    end
end