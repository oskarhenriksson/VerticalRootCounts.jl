using Documenter
using VerticalRootCounts
using Oscar
using Catalyst

DocMeta.setdocmeta!(
    VerticalRootCounts,
    :DocTestSetup,
    quote
        using VerticalRootCounts
        using Oscar
        using Catalyst
    end;
    recursive = true,
)

makedocs(
    sitename = "VerticalRootCounts.jl",
    modules = [VerticalRootCounts],
)

deploydocs(
    repo = "github.com/oskarhenriksson/VerticalRootCounts.jl.git",
)