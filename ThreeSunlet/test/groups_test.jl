using Test
include("groups.jl")

Z3 = FiniteCyclicGroup([3])
Z4 = FiniteCyclicGroup([4])
Z2Z2 = FiniteCyclicGroup([2,2])
Z2Z3Z5 = FiniteCyclicGroup([2,3,5])
test_groups = [Z3, Z4, Z2Z2]

test_identities = [getGroupIdentity(G) for G in test_groups]
true_identities = [[0], [0], [0,0]]
for g in eachindex(test_groups)
    @test test_identities[g] == true_identities[g]
end

getGroupElements(Z2Z2)

g1 = [2]
h1 = getInverse(Z4, g1)
@test groupAdd(Z4, [g1, h1]) == getGroupIdentity(Z4)