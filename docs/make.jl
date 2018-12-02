using Documenter, RationalExpectations

makedocs(sitename = "RationalExpectations",
        authors = "Daniel Coutinho",
        doctest = false,
        pages = Any[
                "Home" => "index.md",
                "Sims" => "sims.md",
                "Klein" => "klein.md",
                "Example" => "example.md"
                ]
)

deploydocs(
    deps = Plots,
    repo = "github.com/danmrc/RationalExpectations.jl.git"
)
