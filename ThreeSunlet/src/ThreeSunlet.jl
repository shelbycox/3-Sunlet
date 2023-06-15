module ThreeSunlet

greet() = print("Hello World! I am the 3-Sunlet package for Julia.")

include("cyclic.jl")
include("groups.jl")
include("hyperplanes.jl")
include("lambda_cyclic.jl")
include("lambda_guess.jl")
include("latexTables.jl")
include("ranks_ineq_gen.jl")
include("row_reduction.jl")
include("sampling.jl")
include("weyl.jl")

end # module ThreeSunlet
