include("../src/toDot.jl")

N = 2^15
Z3 = FiniteCyclicGroup([3])
Z22 = FiniteCyclicGroup([2,2])
Z4 = FiniteCyclicGroup([4])

toDot(Z3, N)

toDot(Z22, N)

toDot(Z4, N)