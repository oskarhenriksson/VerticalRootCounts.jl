export cotransversal_presentation

@doc raw"""
    cotransversal_presentation(A::QQMatrix)

Checks if the rowspan matroid of a rational matrix `A` is cotransversal.
If it is, a cotransversal presentation of the matroid is returned; if it is not, `false` is returned.

"""
function cotransversal_presentation(A::QQMatrix)
    M = matroid_from_matrix_columns(A;check=false).pm_matroid;
    transversality_witness = Polymake.matroid.check_transversality(M)
    if transversality_witness == false
        return false
    else
        return [w.+1 for w in transversality_witness]
    end
end
