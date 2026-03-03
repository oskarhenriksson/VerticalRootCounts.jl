using Catalyst, Oscar
using VerticalRootCounts

# The running example
rn = @reaction_network begin
    k1, X1 --> X2
    k2, X2 --> X1
    k3, 2*X1 + X2 --> 3*X1
end;

# Compute the steady state degree directy from the network
@time sd = steady_state_degree(rn)

# The defining matrices for the steady state system
C, M, L = steady_state_system(rn)

# Generic root count of the steady state system
@time generic_root_count(C, M, L)

# Without the transversality check
@time generic_root_count(C, M, L, check_transversality=false)

# Lower bound directly from the network
@time bound, b, k, h = lower_bound_of_maximal_positive_steady_state_count(rn)
@time bound, b, k, h = lower_bound_of_maximal_positive_root_count(C, M, L)

# Vertify the result for a given choice of b and h
b = [71]
k = [839, 562, 13]
h = [37,97,18]
lower_bound_of_maximal_positive_root_count_fixed_b_k_h(C, M, L, b, k, h)

# Cell cycle
rn = Catalyst.@reaction_network begin
    k1, C + Mp --> C + M
    k2, Cp + M --> C + M
    k3, M + W --> Mp + W
    k4, M + W --> M + Wp
    k5, C --> Cp
    k6, Wp --> W
end;
@time steady_state_degree(rn)
@time lower_bound_of_maximal_positive_steady_state_count(rn)

# Verify the result for a given choice of b and h
C, M, L = steady_state_system(rn)
b =  [69, 42, 81]
k = [622, 732, 905, 631, 567, 253]
h = [12, 86, 11, 27, 84]
lower_bound_of_maximal_positive_root_count_fixed_b_k_h(C, M, L, b, k, h)

# HHK
rn = Catalyst.@reaction_network begin
    k1, HK00 --> HKp0
    k2, HKp0 -->  HK0p
    k3, HK0p --> HKpp  
    k4, HK0p  + Hpt --> HK00 + Hptp
    k5, HKpp  + Hpt --> HKp0 + Hptp
    k6, Hptp  --> Hpt
end;
@time steady_state_degree(rn)
@time lower_bound_of_maximal_positive_steady_state_count(rn)

# Verify the result for a given choice of b and h
C, M, L = steady_state_system(rn)
b = [59, 34]
k = [839, 562, 13, 421, 233, 109]
h = [84, 46, 30, 13, 23, 68]
lower_bound_of_maximal_positive_root_count_fixed_b_k_h(C, M, L, b, k, h)
  

# Triangle network
rn = Catalyst.@reaction_network begin
  k1, 3*X1 + 2*X2 --> 6*X1
  k2, 3*X1 + 2*X2 --> 4*X2 
  k3, 4*X2 --> 3*X1 + 2*X2 
  k4, 6*X1 -->  4*X2
end

C, M, L = steady_state_system(rn)

generic_root_count(C, M, L)
bound, _, _ = lower_bound_of_maximal_positive_root_count(C, M, L)
A = matrix(ZZ, [[3, 2]])
toric_root_bound(A, L)
bound, _, _ = toric_lower_bound_of_maximal_positive_root_count(A, L)


# 1-site phosphorylation
rn = Catalyst.@reaction_network begin
  k1, S0 + E --> ES0
  k2, ES0  --> S0 + E
  k3, ES0  --> S1+E
  k4, S1 + F  --> FS1
  k5, FS1  --> S1 + F
  k6, FS1 --> S0 + F
end;

@time steady_state_degree(rn)
@time lower_bound_of_maximal_positive_steady_state_count(rn)

# Verify the result for a given choice of b and h
C, M, L = steady_state_system(rn)
b = [68, 52, 99]
k = [839, 562, 13, 421, 233, 109]
h = [79, 26, 89, 92]
lower_bound_of_maximal_positive_root_count_fixed_b_k_h(C, M, L, b, k, h)

A = kernel(matrix(ZZ, hcat([M[:, i] - M[:, ncols(M)] for i=1:ncols(M)-1]...)))
toric_root_bound(A, L, check_transversality=true)
toric_lower_bound_of_maximal_positive_root_count(A, L)

#2-site phosphorylation
rn = Catalyst.@reaction_network begin
    @parameters k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12
    @species E(t) F(t)  S0(t) S1(t) ES0(t) FS1(t) S2(t) ES1(t) FS2(t)
  k1, S0 + E --> ES0
  k2, ES0  --> S0 + E
  k3, ES0  --> S1+E
  k4, S1 + F  --> FS1
  k5, FS1  --> S1 + F
  k6, FS1 --> S0 + F
  k7, S1 + E --> ES1
  k8, ES1 --> S1 + E
  k9, ES1 --> S2 + E
  k10, S2 + F  -->FS2
  k11, FS2 --> S2 + F
  k12, FS2 --> S1 + F
end 

@time steady_state_degree(rn)
@time lower_bound_of_maximal_positive_steady_state_count(rn)
b = [1811, 1135, 4769]
k = [744, 59, 746, 120, 270, 517, 377, 798, 632, 431, 722, 333]
h = [259, 800, 750, 684, 363, 120, 524, 616]
C, M, L = steady_state_system(rn)
lower_bound_of_maximal_positive_root_count_fixed_b_k_h(C, M, L, b, k, h)

toric_root_bound(A, L, check_transversality=false)
toric_lower_bound_of_maximal_positive_root_count(A, L)


# Wnt pathway
rn = @reaction_network begin
  k1,  X1          --> X2
  k2,  X2          --> X1
  k3,  X2 + X4     --> X14
  k4,  X14         --> X2 + X4
  k5,  X14         --> X2 + X5
  k6,  X5 + X8     --> X16
  k7,  X16         --> X5 + X8
  k8,  X16         --> X4 + X8
  k9,  X4 + X10    --> X18
  k10, X18         --> X4 + X10
  k11, X18         --> X4
  k12, 0           --> X10
  k13, X10         --> 0
  k14, X3 + X6     --> X15
  k15, X15         --> X3 + X6
  k16, X15         --> X3 + X7
  k17, X7 + X9     --> X17
  k18, X17         --> X7 + X9
  k19, X17         --> X6 + X9
  k20, X6 + X11    --> X19
  k21, X19         --> X6 + X11
  k22, X19         --> X6
  k23, X11         --> 0
  k24, X11 + X12   --> X13
  k25, X13         --> X11 + X12
  k26, X2          --> X3
  k27, X3          --> X2
  k28, X5          --> X7
  k29, X7          --> X5
  k30, X10         --> X11
  k31, X11         --> X10
end;

@time sd = steady_state_degree(rn; verbose=true)

# Example with more isolated zeros than the generic root count
C = matrix(QQ, [
     1  -12   58  -144  193  -132    1  -12   58  -144  193  -132  0   0   0   72  0
     0    0    0     0    0     0    0    0    0     0    0     1  -1   1   0  -1  0
     0    0    0     0    0     0    0    0    0     0    0     0   0   0    2  -2   0
])

M = matrix(ZZ, [
    6  5  4  3  2  1  0  0  0  0  0  0  0  0  0 0  0
    0  0  0  0  0  0  6  5  4  3  2  1  1  0  0 0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  1  1 0  0
])

generic_root_count(C, M)
