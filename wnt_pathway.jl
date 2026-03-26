using Catalyst, Oscar
using VerticalRootCounts

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

F = steady_state_system(rn)
mixed_volume(F)
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
