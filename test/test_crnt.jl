using Test
using Catalyst
using VerticalRootCounts

@testset verbose=true "Steady state systems" begin
    @testset "4-site phosphorylation" begin

        C, M, L, A = multisite_phosphorylation_matrices(4)

        rn = @reaction_network begin
            @species E(t) F(t) S0(t) S1(t) S2(t) S3(t) S4(t) ES0(t) ES1(t) ES2(t) ES3(t) FS1(t) FS2(t) FS3(t) FS4(t)
            k1, S0 + E --> ES0
            k2, ES0 --> S0 + E
            k3, ES0 --> S1 + E
            k4, S1 + F --> FS1
            k5, FS1 --> S1 + F
            k6, FS1 --> S0 + F
            k7, S1 + E --> ES1
            k8, ES1 --> S1 + E
            k9, ES1 --> S2 + E
            k10, S2 + F --> FS2
            k11, FS2 --> S2 + F
            k12, FS2 --> S1 + F
            k13, S2 + E --> ES2
            k14, ES2 --> S2 + E
            k15, ES2 --> S3 + E
            k16, S3 + F --> FS3
            k17, FS3 --> S3 + F
            k18, FS3 --> S2 + F
            k19, S3 + E --> ES3
            k20, ES3 --> S3 + E
            k21, ES3 --> S4 + E
            k22, S4 + F --> FS4
            k23, FS4 --> S4 + F
            k24, FS4 --> S3 + F
        end

        F = steady_state_system(rn)

        @test rref(C)[2] == rref(F.C)[2]
        @test M == F.M
        @test rref(L)[2] == rref(F.L)[2]
    end

end


