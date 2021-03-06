using Documenter, RationalExpectations, Plots

makedocs(sitename = "RationalExpectations",
        authors = "Daniel Coutinho",
        doctest = false,
        pages = Any[
                "Home" => "index.md",
                "Sims" => "sims.md",
                "Klein" => "klein.md",
                "Example" => "example.md",
                "Bibliography"=> "biblio.md"
                ]
)

deploydocs(
    repo = "github.com/danmrc/RationalExpectations.jl.git"
)
