include("hyperplanes.jl")
using Test

Z3 = FiniteCyclicGroup([3])
H3 = [[1, 0, 0, 0, 0], [0, 1, 0, -1, 0], [0, 1, 0, 0, -1], [0, 1, 0, 1, 1], [0, 0, 1, 1, 0], [0, 0, 1, 0, 1], [0, 0, 1, -1, -1]]

@test getMuVector([[0], [0], [0]], Z3) == [1, 0 ,0]
@test getEtaVector([[1], [2], [0]], Z3) == [1, 1]
@test getEtaVector([[0], [0], [0]], Z3) == [0, 0]

@test getHyperplane(([0], [0], [0]), Z3) == [1, 0, 0, 0, 0]
@test getHyperplane(([1], [0], [2]), Z3) == [0, 1, 0, -1, 0]
@test getHyperplane(([1], [1], [1]), Z3) == [0, 1, 0, 0, -1]
@test getHyperplane(([1], [2], [0]), Z3) == [0, 1, 0, 1, 1]

@test Set(generateSunletArr(Z3)) == Set(H3) ## some sign errors here... need to look into this further
print(generateSunletArr(Z3))
print(getValidGroupTriples(Z3))

@test -1*[1,1] == [-1,-1]
@test getEtaSign(0,2) == 1

Z2Z2 = FiniteCyclicGroup([2,2])
generateSunletArr(Z2Z2)