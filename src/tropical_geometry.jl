
@doc raw"""
    perturb_and_intersect_if_transversal(TropL::TropicalLinearSpace, TropB::TropicalVariety; perturbation::Vector{<:Integer}=rand(Int16,ambient_dim(TropL)), with_multiplicities::Bool = true)

Perturb `TropB` by `perturbation` and intersect with `TropL` if the intersection is transversal.  Return four things:
- the vector `perturbation`
- a Boolean `true` or `false` depending on whether `TropB + perturbation` and `TropL` intersect transversally
- a vector that contains the intersection points if the intersection is transversal (otherwise empty)
- a vector that contains the multiplicities of the intersection points if the intersection is transversal and `with_multiplicities` is `true` (otherwise empty)

Assumes that `TropL` is a polyhedral fan and that `TropB` is the tropicalization of a binomial variety, i.e., a regular linear space.
"""
function perturb_and_intersect_if_transversal(TropL::TropicalLinearSpace,
    TropB::TropicalVariety;
    perturbation::Vector{<:Integer}=rand(Int16,ambient_dim(TropL)),
    with_multiplicities::Bool = true
)

    bergmanRays, bergmanLineality = rays_modulo_lineality(TropL)
    bergmanRays = matrix(QQ, bergmanRays)
    bergmanLineality = matrix(QQ, bergmanLineality)

    minimalFaces, linearSpaceBasis = minimal_faces(TropB)
    linearSpaceBasis = matrix(QQ, linearSpaceBasis)

    @req length(minimalFaces) == 1 "Several minimal faces found in TropB"

    # compute the projection matrix onto the orthogonal complement of the euclidean linear space
    basisOfComplementTransposed = kernel(linearSpaceBasis, side=:right)
    basisOfComplement = transpose(basisOfComplementTransposed)
    projectionMatrix = basisOfComplementTransposed * inv(basisOfComplement * basisOfComplementTransposed) * basisOfComplement

    # project the rays of the Bergman fan
    projectedRays = bergmanRays * projectionMatrix
    projectedLineality = bergmanLineality * projectionMatrix

    # make it consistent whether projectionPerturbation and perturbation are rows/columns
    projectedPerturbation = matrix(QQ, [perturbation]) * projectionMatrix
    stableIntersectionPoints = Vector{QQFieldElem}[]
    stableIntersectionMults = Int[]

    indicesOfCones = ray_indices(maximal_polyhedra(TropL))
    nRaysPerCone = sum(indicesOfCones[1, :])
    for i in 1:nrows(indicesOfCones)
        # read off rays of the projected cone
        indicesOfCone = findall(indicesOfCones[i, :])
        projectedRaysOfCone = projectedRays[indicesOfCone, :]

        # test whether projected direction lies in projected cone
        # warning: be careful about the sign of the perturbation
        can_solve, solution = can_solve_with_solution(vcat(projectedRaysOfCone, projectedLineality),
            projectedPerturbation; side=:left)
        if can_solve
            firstZero = findfirst(isequal(0), solution)
            if (firstZero != nothing) && (firstZero[2] <= nRaysPerCone)
                # random direction not generic, lies on the boundary of the cone
                return (
                    perturbation = perturbation,
                    is_transversal = false,
                    points = stableIntersectionPoints,
                    multiplicities = stableIntersectionMults
                )
            end
            firstNegative = findfirst(a -> (a < 0), solution)
            if (firstNegative == nothing) || (firstNegative[2] > nRaysPerCone)
                # random direction lies in the interior of the cone,
                # compute intersection point and intersection multiplicity
                intersectionPoint = solution * vcat(bergmanRays[indicesOfCone, :], bergmanLineality)

                push!(stableIntersectionPoints, intersectionPoint[1, :])
                coneSpanBasis = vcat(bergmanRays[indicesOfCone, :], bergmanLineality)
                if with_multiplicities == true
                    push!(stableIntersectionMults, tropical_intersection_multiplicity(coneSpanBasis, linearSpaceBasis))
                end
            end
        end
    end
 
    return (
        perturbation = perturbation,
        is_transversal = true,
        points = stableIntersectionPoints,
        multiplicities = stableIntersectionMults
    )
end


function tropical_intersection_multiplicity(B1, B2)
    @assert ncols(B1) == ncols(B2) && nrows(B1) + nrows(B2) >= ncols(B1)

    # primitive scales every row by the lcm of the denominators, making the matrix integral
    # saturate computes a basis of the saturation of the sublattice spanned by the row vectors
    B1 = saturate(matrix(ZZ, Polymake.common.primitive(B1)))
    B2 = saturate(matrix(ZZ, Polymake.common.primitive(B2)))

    return abs(prod(elementary_divisors(vcat(B1,B2)))) |> Int
end
