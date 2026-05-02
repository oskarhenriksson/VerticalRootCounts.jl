export steady_state_system,
    multisite_phosphorylation_matrices

"""
    steady_state_system(rn)

    For a reaction network `rn`, compute the coefficient, exponent, and linear part matrices
    of the augmented vertically parametrized steady state system.

    # Example
    ```jldoctest
    julia> using Catalyst;
    
    julia> rn = @reaction_network begin
        k1, X1 --> X2
        k2, X2 --> X1
        k3, 2*X1 + X2 --> 3*X1
    end;

    julia> F = steady_state_system(rn);

    julia> F.C
    [1   -1   -1]

    julia> F.M
    [1   0   2]
    [0   1   1]

    julia> F.L
    [1   1]
    ```

"""
function steady_state_system(rn::ReactionSystem)

    @req all(ismassaction.(reactions(rn), Ref(rn))) "All reactions must have mass action kinetics"

    N = matrix(QQ, netstoichmat(rn)) #stoiciometric matrix
    M = matrix(ZZ, substoichmat(rn)) #kinetic matrix
    L = rref(matrix(QQ, conservationlaws(rn)))[2] #conserved quantities
    
    # Compute the matrix C (choose of linearly independent rows of N)
    first_nonzero_indices = [findfirst(!iszero, row) for row in eachrow(L)]
    C = N[setdiff(collect(1:nrows(N)), first_nonzero_indices), :] 

    return AugmentedVerticalSystem(C, M, L)
end





function multisite_phosphorylation_matrices(k)

    # Ordering of species: E, F, S0, ..., Sm, ES0, ..., ESm-1, FS1, ..., FSm

    # Stoichiometric matrix
    e_row = vcat([[-1, 1, 1, 0, 0, 0] for i in 1:k]...) |> transpose
    f_row = vcat([[0, 0, 0, -1, 1, 1] for i in 1:k]...) |> transpose
    S_rows = zeros(Int, k + 1, 6 * k)
    for i = 1:k
        S_rows[i:i+1, 6*(i-1)+1:6*i] = [-1 1 0 0 0 1; 0 0 1 -1 1 0]
    end
    ES_rows = zeros(Int, k, 6 * k)
    for i = 1:k
        ES_rows[i, 6*(i-1)+1:6*i] = [1, -1, -1, 0, 0, 0]
    end
    FS_rows = zeros(Int, k, 6 * k)
    for i = 1:k
        FS_rows[i, 6*(i-1)+1:6*i] = [0, 0, 0, 1, -1, -1]
    end
    N = matrix(QQ, vcat(e_row, f_row, S_rows, ES_rows, FS_rows))

    # Kinetic matrix
    e_rows = vcat([[1, 0, 0, 0, 0, 0] for i in 1:k]...) |> transpose
    f_rows = vcat([[0, 0, 0, 1, 0, 0] for i in 1:k]...) |> transpose
    S_rows = zeros(Int, k + 1, 6 * k)
    for i = 1:k
        S_rows[i:i+1, 6*(i-1)+1:6*i] = [1 0 0 0 0 0; 0 0 0 1 0 0]
    end
    ES_rows = zeros(Int, k, 6 * k)
    for i = 1:k
        ES_rows[i, 6*(i-1)+1:6*i] = [0, 1, 1, 0, 0, 0]
    end
    FS_rows = zeros(Int, k, 6 * k)
    for i = 1:k
        FS_rows[i, 6*(i-1)+1:6*i] = [0, 0, 0, 0, 1, 1]
    end
    M = matrix(ZZ, vcat(e_rows, f_rows, S_rows, ES_rows, FS_rows))

    # Conservation laws
    s_tot = [0; 0; ones(Int, 3 * k + 1)] |> transpose
    e_tot = [1; 0; zeros(Int, k+1); ones(Int, k); zeros(Int, k)] |> transpose
    f_tot = [0; 1; zeros(Int, k+1); zeros(Int, k); ones(Int, k)] |> transpose
    L = matrix(QQ, vcat(s_tot, e_tot, f_tot))

    @assert all(iszero, L*N) "Conservation laws matrix L is not orthogonal to stoichiometric matrix N" 

    # Coefficient matrix
    first_nonzero_indices = [findfirst(!iszero, row) for row in eachrow(L)]
    C = N[setdiff(collect(1:nrows(N)), first_nonzero_indices), :] 

    # Exponent matrix for parametric toricity
    row1 = [1; 0; 0:k; 1:k; 1:k] |> transpose
    row2 = [0; 1; -(0:k); -(0:k-1); -(0:k-1)] |> transpose
    row3 = [0; 0; ones(Int, 3*k+1)] |> transpose
    A = matrix(ZZ, vcat(row1, row2, row3))

    return C, M, L, A
end