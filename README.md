# VerticalRootCounts.jl
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://oskarhenriksson.github.io/VerticalRootCounts.jl/dev/)
[![CI](https://github.com/oskarhenriksson/VerticalRootCounts.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/oskarhenriksson/VerticalRootCounts.jl/actions/workflows/ci.yml)

This Julia package is a proof-of-concept implementation of the tropical root bounds presented in the forthcoming preprint **Root bounds of vertical system using tropical geometry** by Elisenda Feliu, Paul Helminck, Oskar Henriksson, Yue Ren, Benjamin Schröter, and Máté L. Telek.

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
 Computation method: mixed volume for cotransversal presentation
 Choice of parameters a: [447, 884, 710, 139, 73, 644]
 Choice of constant terms b: [2132, 438, 369]
 Row supports of cotransversal presentation for the nonlinear part: 
  [1, 2, 3, 4]
  [1, 2, 3, 4]
  [1, 2, 3, 4]
 Row supports of cotransversal presentation for the linear part: 
  [1, 2, 3, 4, 5, 6, 7]
  [5, 6, 7]
  [2, 3, 7]

julia> generic_root_count(F; check_cotransversality=false)
Result of generic root count computation
========================================
 Generic root count: 3
 Computation method: stable intersection of binomial and linear parts
 Choice of parameters a: [298, 748, 393, 553, 853, 358]
 Choice of constant terms b: [2715, 1609, 1198]
 Choice of perturbation h: [20049, 14768, 16258, 25025, -24345, 5994, -14831, -28978, 17231, 31370, 26329]

```

We can also compute a lower bound on the maximal numer of positive zeros. 

```julia-repl
julia> lower_bound_of_maximal_positive_root_count(F)
Result of positive tropical root bound computation
==================================================
 Lower bound on the maximal number of positive roots: 1
 Computation method: stable intersection of binomial and linear parts
 Choice of parameters a: [572, 117, 742, 551, 18, 846]
 Choice of constant terms b: [1750, 1608, 1024]
 Choice of perturbation h: [873, 352, 327, 768]

```

To make use of toricity with respect to a known exponent matrix `A`, we instead use the `toric_root_bound` command:

```julia-repl
julia> A = matrix(ZZ, [1 0 1 0 1 1; 0 1 1 0 1 1; 0 0 0 1 -1 0]);

julia> toric_root_bound(A, F)
Result of toric root bound computation
======================================
 Toric root bound: 3
 Computation method: mixed volume for cotransversal presentation
 Choice of constant terms b: [47751, 28927, -11192]
 Row supports of cotransversal presentation for the linear part: 
  [1, 2, 3, 4, 5, 6, 7]
  [5, 6, 7]
  [2, 3, 7]


julia> toric_lower_bound_of_maximal_positive_root_count(A, F)
Result of positive toric root bound computation
===============================================
 Lower bound on the maximal number of positive roots: 1
 Computation method: stable intersection of binomial and linear parts
 Choice of constant terms b: [2688, 791, 1540]
 Choice of perturbation h: [700, 288, 806, 425, 969, 791, 939]

```

We can also apply techniques directly to chemical reaction networks through our interface to `Catalyst.jl` through the following commands:

```julia-repl
julia> using Catalyst;
julia> rn = Catalyst.@reaction_network begin
         k1, S0 + K --> KS0
         k2, KS0  --> S0 + K
         k3, KS0  --> S1+K
         k4, S1 + P  --> PS1
         k5, PS1  --> S1 + P
         k6, PS1 --> S0 + P
       end;

julia> steady_state_degree(rn)
Result of generic root count computation
========================================
 Generic root count: 3
 Computation method: mixed volume for cotransversal presentation
 Choice of parameters a: [397, 604, 976, 64, 93, 94]
 Choice of constant terms b: [2290, 627, 1135]
 Row supports of cotransversal presentation for the nonlinear part: 
  [1, 2, 3, 4]
  [1, 2, 3, 4]
  [1, 2, 3, 4]
 Row supports of cotransversal presentation for the linear part: 
  [1, 2, 3, 4, 5, 6, 7]
  [5, 6, 7]
  [2, 3, 7]

julia> lower_bound_of_maximal_positive_steady_state_count(rn)
Result of positive tropical root bound computation
==================================================
 Lower bound on the maximal number of positive roots: 1
 Computation method: stable intersection of binomial and linear parts
 Choice of parameters a: [51, 918, 595, 982, 625, 346]
 Choice of constant terms b: [1762, 534, 534]
 Choice of perturbation h: [366, 602, 750, 30]

```
