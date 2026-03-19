using Test
using Random, Oscar, Catalyst
using VerticalRootCounts

@testset verbose=true "Root bounds for networks" begin
@testset "Small example" begin

    Random.seed!(1234)

    C = matrix(QQ, [[1,-1,-1]])
    M = matrix(ZZ, [[1,0,2], [0,1,1]])
    L = matrix(QQ, [[1,1]])

    F = AugmentedVerticalSystem(C, M, L)

    @test generic_root_count(F).count == 3
    @test has_nondegenerate_zero(F)
    
    rn = @reaction_network begin
        k1, X1 --> X2
        k2, X2 --> X1
        k3, 2*X1 + X2 --> 3*X1
    end;
    
    @test steady_state_degree(rn, check_cotransversality=true).count == 3
    @test steady_state_degree(rn, check_cotransversality=false).count == 3

    a = [83, 56, 13]
    b = [71]
    h = [37,97,18]
    @test lower_bound_of_maximal_positive_root_count_fixed_a_b_h(F, a, b, h).bound == 3

    @test lower_bound_of_maximal_positive_root_count(F, num_a_b_attempts=5, num_h_attempts_per_a_b=5).bound == 3
    
    @test lower_bound_of_maximal_positive_steady_state_count(rn, num_a_b_attempts=3, num_h_attempts_per_a_b=3).bound == 3
    
end

@testset "Degenerate purely vertical system" begin
    
    Random.seed!(1234)

    C = matrix(QQ, [-1  0  0  0  1  0; 0 -1  0  0  0  1;  0  0 -1  1  0  0]);
    M = matrix(ZZ, [3  2  1  0  0  0; 0  1  0  2  1  0; 0  0  1  0  1  2]);

    F = AugmentedVerticalSystem(C, M)
    
    @test !has_nondegenerate_zero(F)
    
    @test generic_root_count(F, check_cotransversality=true).count == 0
    @test generic_root_count(F, check_cotransversality=false).count == 0
    @test_nowarn sprint(show, MIME("text/plain"), generic_root_count(F, check_cotransversality=true))
    @test lower_bound_of_maximal_positive_root_count(F).bound == 0
    @test_nowarn sprint(show, MIME("text/plain"), lower_bound_of_maximal_positive_root_count(F))


end

@testset "Cell cycle" begin

    Random.seed!(1234)

    rn = @reaction_network begin
        k1, C + Mp --> C + M
        k2, Cp + M --> C + M
        k3, M + W --> Mp + W
        k4, M + W --> M + Wp
        k5, C --> Cp
        k6, Wp --> W
    end

    @test steady_state_degree(rn).count == 2

    F = steady_state_system(rn)
    a = QQ.(1//3 * [284, 215, 921, 770, 883, 792])
    b = QQ.(1//5 * [69, 42, 81])
    h = [12, 86, 11, 27, 84]
    @test lower_bound_of_maximal_positive_root_count_fixed_a_b_h(F, a, b, h).bound == 2

end

@testset "HHK network" begin

    Random.seed!(1234)
    
    rn = @reaction_network begin
        k1, HK00 --> HKp0
        k2, HKp0 -->  HK0p
        k3, HK0p --> HKpp  
        k4, HK0p  + Hpt --> HK00 + Hptp
        k5, HKpp  + Hpt --> HKp0 + Hptp
        k6, Hptp  --> Hpt
    end

    @test steady_state_degree(rn).count == 3

    F = steady_state_system(rn)
    a = [84, 46, 30, 13, 23, 68]
    b = [59, 34]
    h = [834, 131, 91, 217, 253, 498]
    @test lower_bound_of_maximal_positive_root_count_fixed_a_b_h(F, a, b, h).bound == 3

end

@testset "1-site phosphorylation" begin

    Random.seed!(1234)

    rn = @reaction_network begin
        k1, S0 + E --> ES0
        k2, ES0  --> S0 + E
        k3, ES0  --> S1+E
        k4, S1 + F  --> FS1
        k5, FS1  --> S1 + F
        k6, FS1 --> S0 + F
    end

    @test steady_state_degree(rn).count == 3

    F = steady_state_system(rn)

    @test generic_degree(AugmentedVerticalSystem(F.C, F.M)).count == 4

    a = [84, 46, 30, 13, 23, 68]
    b = [68, 52, 99]
    h = [79, 26, 89, 92]
    @test lower_bound_of_maximal_positive_root_count_fixed_a_b_h(F, a, b, h).bound == 1

    M = F.M
    A = kernel(matrix(ZZ, hcat([M[:, i] - M[:, ncols(M)] for i=1:ncols(M)-1]...)))
    @test toric_root_bound(A, F, check_cotransversality=true).bound == 3
    @test toric_root_bound(A, F, check_cotransversality=false).bound == 3

    h = [936, 145, 170, 323, 169, 271, 439]
    @test toric_lower_bound_of_maximal_positive_root_count_fixed_b_h(A, F, b, h).bound == 1

end

@testset "2-site phosphorylation" begin

    Random.seed!(1234)

    rn = @reaction_network begin
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

    @test steady_state_degree(rn).count == 5

end


@testset "Triangle network" begin

    Oscar.set_seed!(1234)

    rn = Catalyst.@reaction_network begin
        k1, 3*X1 + 2*X2 --> 6*X1
        k2, 3*X1 + 2*X2 --> 4*X2 
        k3, 4*X2 --> 3*X1 + 2*X2 
        k4, 6*X1 -->  4*X2
    end;

    F = steady_state_system(rn)

    @test generic_root_count(F).count == 6
    @test lower_bound_of_maximal_positive_root_count(F).bound == 1
    A = matrix(ZZ, [[3, 2]])
    @test toric_root_bound(A, F).bound == 3
    @test toric_lower_bound_of_maximal_positive_root_count(A, F).bound == 1

end

@testset "Mixed volume for Wnt system" begin
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
    @test mixed_volume(F) == 56

end

end
