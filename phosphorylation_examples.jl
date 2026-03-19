using VerticalRootCounts

for k = 1:5
    println("\n\n$(k)-site phosphorylation network:\n")
    C, M, L, A = multisite_phosphorylation_matrices(k)
    F = AugmentedVerticalSystem(C, M, L)
    if k<=2
        t_grc = @elapsed grc = generic_root_count(F, check_transversality=false).count
        println("Generic root count: $grc (computed in $t_grc seconds)")
        t_pos = @elapsed positive_bound = lower_bound_of_maximal_positive_root_count(F; show_progress=false).bound
        println("Lower bound of maximal positive steady state count: $positive_bound (computed in $t_pos seconds)")
    else
        println("Generic root count: skipped")
        println("Lower bound of maximal positive steady state count: skipped")
    end
    t_toric = @elapsed grc_toric = toric_root_bound(A, F, check_transversality=false).bound
    println("Toric steady state bound: $grc_toric (computed in $t_toric seconds)")
    t_toric_pos = @elapsed toric_positive_bound = toric_lower_bound_of_maximal_positive_root_count(A, F).bound
    println("Toric lower bound of maximal positive steady state count: $toric_positive_bound (computed in $t_toric_pos seconds)")
    t_grc_mat = @elapsed grc_mat = generic_root_count(F, check_transversality=true).count
    println("Generic root count with transversality: $grc_mat (computed in $t_grc_mat seconds)")
    t_toric_mat = @elapsed grc_toric_mat = toric_root_bound(A, F, check_transversality=true).bound
    println("Toric steady state bound with transversality: $grc_toric_mat (computed in $t_toric_mat seconds)")
end