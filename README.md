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

For the running example in the paper, we compute the generic root count as follows:

```julia-repl
julia> using Oscar;

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
 Choice of constant terms b: [1733, 824, 1175]
 Choice of parameters k: [154, 882, 554, 467, 512, 256]
 Computation method: mixed volume

```

The default is to try to compute the generic root count as a mixed volume if the cotransversality result applies. If it does not apply (or if we turn it off through the `check_cotransversality` keyword), we instead use the main theorem.

```julia-repl
julia> generic_root_count(F, check_cotransversality=false)
Result of generic root count computation
========================================
 Generic root count: 3
 Choice of constant terms b: [2127, 1106, 1075]
 Choice of parameters k: [19, 992, 670, 353, 833, 424]
 Computed method: stable intersection
 Computation of perturbation h: [-31132, 27262, -23936, 25428, -27814, 31751, -19714, -28295, -16342, -25071, 8234]

```

We can also compute a lower bound on the maximal numer of positive zeros. 

```julia-repl
julia> lower_bound_of_maximal_positive_root_count(F)
Result of positive tropical root bound computation
==================================================
 Lower bound on the maximal number of positive roots: 1
 Choice of constant terms b: QQFieldElem[2152, 885, 1104]
 Choice of parameters k: [370, 833, 479, 588, 369, 149]
 Choice of perturbation h: [359, 87, 248, 609]
```

To make use of toricity with respect to a known exponent matrix `A`, we instead use the `toric_root_bound` command:

```julia-repl
julia> A = matrix(ZZ, [1 0 1 0 1 1; 0 1 1 0 1 1; 0 0 0 1 -1 0]);

julia> toric_root_bound(A, F)
Result of toric root bound computation
======================================
 Toric root bound: 3
 Choice of constant terms b: [85081, 35610, 17658]
 Computation method: mixed volume
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
 Choice of constant terms b: [3268, 1638, 1665]
 Choice of parameters k: [506, 26, 187, 315, 829, 561]
 Computation method: mixed volume

julia> lower_bound_of_maximal_positive_steady_state_count(rn)
Result of positive tropical root bound computation
==================================================
 Lower bound on the maximal number of positive roots: 1
 Choice of constant terms b: QQFieldElem[2958, 1490, 664]
 Choice of parameters k: [137, 676, 37, 172, 671, 273]
 Choice of perturbation h: [435, 591, 298, 691]


```

For further examples, we refer to the Jupyter notebook `examples.ipynb`.
