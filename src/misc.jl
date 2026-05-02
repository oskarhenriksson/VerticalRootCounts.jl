export standard_vector,
mixed_volume

standard_vector = (i, n) -> [j == i ? 1 : 0 for j in 1:n]


"""
    mixed_volume(F::Vector{<:MPolyRingElem})

Computes a the mixed volume of a square polynomial system (encoded as a vector of polynomials).

"""
function mixed_volume(F::Vector{<:MPolyRingElem})
    @assert all(parent(f) == parent(first(F)) for f in F) "All polynomials need to be in the same polynomial ring"
    @assert length(F) == nvars(parent(first(F))) "The number of polynomials needs to be equal to the number of variables"
    supports = Matrix{Int}[hcat(exponents(f)...) for f in F]
    return mixed_volume(supports)
end


mixed_volume(F::AugmentedVerticalSystem) = mixed_volume(F.system)