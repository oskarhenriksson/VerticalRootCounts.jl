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

For the running example in the paper, we compute the generic root count as follows.
The default is to try to compute the generic root count as a mixed volume if the cotransversality result applies. If it does not apply (or if we turn it off through the `check_cotransversality` keyword), we instead use the main theorem.

```julia-repl
julia> using Oscar, VerticalRootCounts

julia> C = matrix(QQ, [1 -1 -1 0 0 0; 0 0 1 -1 1 0; 0 0 0 1 -1 -1]);

julia> M = matrix(ZZ, [1 0 0 0 0 0; 1 0 0 0 0 0; 0 1 1 0 0 0; 0 0 0 1 0 0; 0 0 0 1 0 0; 0 0 0 0 1 1]);

julia> L = matrix(QQ, [1 0 1 1 0 1; 0 1 1 0 0 0; 0 0 0 0 1 1]);

julia> F = AugmentedVerticalSystem(C, M, L)
Augmented vertical system
=========================
Variables: x[1], x[2], x[3], x[4], x[5], x[6]
Parameters: a[1], a[2], a[3], a[4], a[5], a[6], b[1], b[2], b[3]

 a[1]*x[1]*x[2] + (-a[2] - a[3])*x[3]
 a[3]*x[3] - a[4]*x[4]*x[5] + a[5]*x[6]
 a[4]*x[4]*x[5] + (-a[5] - a[6])*x[6]
 x[1] + x[3] + x[4] + x[6] - b[1]
 x[2] + x[3] - b[2]
 x[5] + x[6] - b[3]

julia> generic_root_count(F)
Result of generic root count computation
========================================
 Generic root count: 3
 Choice of parameters a: [982, 332, 647, 886, 866, 326]
 Choice of constant terms b: [2110, 1837, 826]
 Computation method: mixed volume of cotransversal presentation

julia> generic_root_count(F; check_cotransversality=false)
Result of generic root count computation
========================================
 Generic root count: 3
 Choice of parameters a: [836, 343, 970, 876, 458, 272]
 Choice of constant terms b: [1667, 826, 1664]
 Computation method: stable intersection of binomial and linear parts
 Computation of perturbation h: [-17615, -18571, -12785, -9616, -23690, -12039, -914, 3718, 615, 16508, 30700]

```

We can also compute a lower bound on the maximal numer of positive zeros. 

```julia-repl
julia> lower_bound_of_maximal_positive_root_count(F)
Result of positive tropical root bound computation
==================================================
 Lower bound on the maximal number of positive roots: 1
 Choice of parameters a: [271, 779, 555, 109, 770, 460]
 Choice of constant terms b: [1319, 1004, 837]
 Choice of perturbation h: [271, 779, 898, 865]
```

To make use of toricity with respect to a known exponent matrix `A`, we instead use the `toric_root_bound` command:

```julia-repl
julia> A = matrix(ZZ, [1 0 1 0 1 1; 0 1 1 0 1 1; 0 0 0 1 -1 0]);

julia> toric_root_bound(A, F)
Result of toric root bound computation
======================================
 Toric root bound: 3
 Choice of constant terms b: [18844, -47913, -6635]
 Computation method: mixed volume for cotransversal presentation

julia> julia> toric_root_bound(A, F; check_cotransversality=false)
Result of toric root bound computation
======================================
 Toric root bound: 3
 Choice of constant terms b: [-16856, -2884, 19205]
 Computation method: stable intersection of binomial and linear parts
 Choice of perturbation h: [-28967, 16168, 20140, 28439, 7112, 10476, -893]

```

We can also apply techniques directly to chemical reaction networks through our interface to `Catalyst.jl` through the following commands:

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
Result of generic root count computation
========================================
 Generic root count: 3
 Choice of parameters a: [854, 96, 139, 19, 404, 89]
 Choice of constant terms b: [2584, 951, 1585]
 Computaion method: mixed volume of cotransversal presentation

julia> lower_bound_of_maximal_positive_steady_state_count(rn)
Result of positive tropical root bound computation
==================================================
 Lower bound on the maximal number of positive roots: 1
 Choice of parameters a: [724, 851, 433, 573, 189, 272]
 Choice of constant terms b: [1211, 1284, 828]
 Choice of perturbation h: [463, 297, 715, 564]


```

For further examples, we refer to the Jupyter notebook `examples.ipynb`.
