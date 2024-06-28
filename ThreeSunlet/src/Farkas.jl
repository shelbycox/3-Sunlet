using Oscar
include("hyperplanes.jl")

G = FiniteCyclicGroup([6])
H = permutedims(hcat(generateSunletArr(G)...))
n = length(H[:,1])
d = length(H[1,:])

signsList = collect(Iterators.flatten(Iterators.product([[1, -1] for i=1:n]...)));
signVecs = [signsList[i:i+n-1] for i=1:n:length(signsList)]

admissibleSigns = []
count = 0
for s in signVecs
    if count % 1000 == 0
        print('.')
    end
    sH = diagm(s)*H
    P = polyhedron((sH, zeros(n)))
    if dim(P) == d
        push!(admissibleSigns, s)
    end
    count = count + 1
end
admissibleSigns

## check that we get the same number
## if same, then we're done
## if different, actually find the ranks of A_lambda