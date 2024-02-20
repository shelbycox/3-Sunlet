include("../src/weyl.jl")
include("../src/lambda_cyclic.jl")

n = 3

v = p(n)
H = getH(n)
O = getOrbit(v, H)
sample = getProj(O, n)
ranks = [getRank(s[1:n], s[(n+1):(2*n)], n) for s in sample]

print(ranks[2], '\n')
for i=eachindex(ranks)
    if isNeighbor(sample[2], sample[i], n)
        print(ranks[i])
    end
end