module ThreeSunlet

greet() = print("Hello World! I am the 3-Sunlet package for Julia.")

## files that need no include statement
include("groups.jl")
include("cyclic.jl")
include("weyl.jl")

## files that need "groups.jl"
include("hyperplanes.jl")
include("lambda_cyclic.jl")
include("lambda_guess.jl")

## files that need some include statements
include("sampling.jl")
include("latexTables.jl")
include("ranks_ineq_gen.jl")

## needs a lot of includes
include("row_reduction.jl")

end # module ThreeSunlet
