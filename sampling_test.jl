using Test
include("sampling.jl")
include("groups.jl")

Z3 = FiniteCyclicGroup([3])

p = getSamplePoint(Z3)
@test length(p) == 5