include("../src/sampling.jl")
include("../src/lambda_guess.jl")
include("../src/lambda_cyclic.jl")

l = 5

G = FiniteCyclicGroup([l+1]) 
H = generateSunletArr(G) 
S = sample(G, H, 2^30) ## (at most) one sample per chamber
C = collect(keys(S)) ## chambers of the arrangement (without duplicates)

maxRank = 5*l + 1

ranks = Dict()
for I in C
   ranks[I] = getRank(S[I], G)
end


global locallyMaximalChamberCount = 0

for I in C
	maximal = true
	rank = ranks[I]
	if rank < maxRank
		for J in C
			if areNeighbors(I, J) && ranks[J] > rank
				maximal = false
			end
		end
		if maximal
            println(I, "\t", rank)
            global locallyMaximalChamberCount += 1
		end
	end
end
println(locallyMaximalChamberCount, "/", length(C), " locally maximal chambers")
