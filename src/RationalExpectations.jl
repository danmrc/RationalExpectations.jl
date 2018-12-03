module RationalExpectations

    using LinearAlgebra

    include("aux.jl")
    include("sims.jl")
    include("klein.jl")

    export
        sims, klein, irf

end #module
