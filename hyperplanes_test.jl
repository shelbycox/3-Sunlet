include("hyperplanes.jl")
using Test

## so far this only works for groups with a single factor

@test getMuVector([[0], [0], [0]], Z3) == [1, 0 ,0]
@test getEtaVector([[1], [2], [0]], Z3) == [-1, -1]
@test getEtaVector([[0], [0], [0]], Z3) == [0, 0]

Z3 = FiniteCyclicGroup([3])
H3 = [[1, 0, 0, 0, 0], [0, 1, 0, -1, 0], [0, 1, 0, 0, -1], [0, 1, 0, -1, -1], [0, 0, 1, 1, 0], [0, 0, 1, 0, 1], [0, 0, 1, 1, 1]]

@test Set(generateSunletArr(Z3)) == Set(H3) ## some sign errors here... need to look into this further

Set(generateSunletArr(Z3))
for T in getValidGroupTriples(Z3)
    println(T)
end