include("groups.jl")
include("lambda_cyclic.jl")
using LinearAlgebra

Z3 = FiniteCyclicGroup([3])
test_lambda = Dict([(4,[0])=>1, (4,[1])=>2, (4,[2])=>-1, (5,[0])=>-2, (5,[1])=>1, (5,[2])=>0, (6,[0])=>4, (6,[1])=>-1, (6,[2])=>2])
M = getMatrix(test_lambda, Z3)
@test size(M) == (15, 9)
@test LinearAlgebra.rank(M) == 9

mu_eta = [0.31276590345517175, 0.4749789656122503, 0.38335449751339157, 0.0032489686671555162, -0.37964053839039946]
test_el = etaToLambda(mu_eta[4:5], getGroupElements(Z3))
@test length(test_el) == 3
parochial_lambda = to_lambda(mu_eta[1:3], mu_eta[4:5], Z3)
@test length(parochial_lambda) == 9
