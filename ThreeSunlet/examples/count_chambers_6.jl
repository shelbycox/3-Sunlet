# using Pkg
# Pkg.activate("ThreeSunlet")
using Oscar
include("../src/hyperplanes.jl")

G = FiniteCyclicGroup([6])
H = permutedims(hcat(generateSunletArr(G)...))
n = length(H[:,1])
d = length(H[1,:])

signsList = collect(Iterators.flatten(Iterators.product([[1, -1] for i=1:n]...)));
signVecs = [signsList[i:i+n-1] for i=1:n:length(signsList)]

admissibleSigns = []
for i=eachindex(signVecs)
    s = signVecs[i]
    if i % 100000 == 0
        print('.')
    end
    sH = diagm(s)*H
    P = polyhedron((sH, zeros(n)))
    if dim(P) == d
        push!(admissibleSigns, s)
    end
end

println()
print(length(admissibleSigns))