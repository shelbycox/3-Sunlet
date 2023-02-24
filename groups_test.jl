using Test
include("groups.jl")

Z3 = FiniteCyclicGroup([3])
Z4 = FiniteCyclicGroup([4])
Z2Z2 = FiniteCyclicGroup([2,2])
test_groups = [Z3, Z4, Z2Z2]

test_identities = [getGroupIdentity(G) for G in test_groups]
true_identities = [[0], [0], [0,0]]
for g in eachindex(test_groups)
    @test test_identities[g] == true_identities[g]
end