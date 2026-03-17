# VerticalRootCounts.jl

[![CI](https://github.com/oskarhenriksson/VerticalRootCounts.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/oskarhenriksson/VerticalRootCounts.jl/actions/workflows/ci.yml)

This Julia package is a proof-of-concept implementation of the tropical root bounds presented in the forthcoming preprint **Tropical bounds of vertical systems** by Elisenda Feliu, Paul Helminck, Oskar Henriksson, Yue Ren, Benjamin Schröter, and Máté L. Telek.

The package is based on tropical functionality from the computer algebra package [Oscar.jl](https://github.com/oscar-system/Oscar.jl/) and the tropical homotopy continuation algorithm for mixed volumes from [MixedSubdivisions.jl](https://github.com/saschatimme/MixedSubdivisions.jl). The package also has an interface to the chemical reaction networks theory package [Catalyst.jl](https://github.com/SciML/Catalyst.jl).

## Installation
To install the package, run the following commands:

```julia
using Pkg
Pkg.add(url="https://github.com/oskarhenriksson/VerticalRootCounts.jl")
```

## Examples of usage

You load the package in a Julia session by running the following command:
```julia
using VerticalRootCounts
```

You can either compute root bounds directly for an augmented vertically parametrized system (defined by a coefficient matrix `C`, an exponent matrix `M`, and coefficient matrix of linear forms `L`) in Oscar format, or for a chemical reaction network given in Catalyst format.

For the 1-site phosphorylation network (Example 2.2 in the paper), we get the steady state degree as follows:

```julia-repl
julia> using Catalyst;
julia> rn = Catalyst.@reaction_network begin
         k1, S0 + E --> ES0
         k2, ES0  --> S0 + E
         k3, ES0  --> S1+E
         k4, S1 + F  --> FS1
         k5, FS1  --> S1 + F
         k6, FS1 --> S0 + F
       end;

julia> steady_state_degree(rn)
3
```

To get a lower bound on the number of positive steady states, we instead run this:

```julia-repl
julia> lower_bound_of_maximal_positive_steady_state_count(rn)
(1, QQFieldElem[901, 977, 1201], [970, 395, 93, 135, 418, 886], [794, 547, 6, 821])

```

The first output is the bound. The remaining outputs are choices of parameters that certify the bound. More precisely, it a choice of total amounts `b`, rate constants `k` and shift vector `h` that give a tropical intersecction with the bound many positve points. (The exact values of these choices depend on ranomized choices inside the algorithm and will therefore vary between runs.)

If we instead want to work directly with the defining matrices of an augmented vertically parametrized system, we can run the following commands:

```julia-repl
julia> using Oscar;

julia> C = matrix(QQ, [1 -1 -1 0 0 0; 0 0 1 -1 1 0; 0 0 0 1 -1 -1]);

julia> M = matrix(ZZ, [1 0 0 0 0 0; 1 0 0 0 0 0; 0 1 1 0 0 0; 0 0 0 1 0 0; 0 0 0 1 0 0; 0 0 0 0 1 1]);

julia> L = matrix(QQ, [1 0 1 1 0 1; 0 1 1 0 0 0; 0 0 0 0 1 1]);

julia> generic_root_count(C, M, L)
3

julia> lower_bound_of_maximal_positive_root_count(C, M, L)
(1, QQFieldElem[1449, 1132, 1538], [696, 838, 259, 713, 11, 312], [142, 958, 851, 938])
```

To make use of toricity with respect to a known exponent matrix `A`, we instead use the `toric_root_bound` command:

```julia-repl
julia> A = matrix(ZZ, [1 0 1 0 1 1; 0 1 1 0 1 1; 0 0 0 1 -1 0]);

julia> toric_root_bound(A, L)
3
```

For further examples, we refer to the Jupyter notebook `examples.ipynb`.
