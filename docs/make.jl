using Documenter, RationalExpectations

makedocs(sitename = "RationalExpectations",
        authors = "Daniel Coutinho",
        doctests = false,
        pages = Any[
                "Home" => "index.md",
                "Sims" => "sims.md",
                "Klein" => "klein.md",
                "Example" => "example.md"
                ]
)
