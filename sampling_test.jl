using Test
include("sampling.jl")
include("groups.jl")
include("hyperplanes.jl")

Z3 = FiniteCyclicGroup([3])
H3 = generateSunletArr(Z3)

p3 = getSamplePoint(Z3)
@test length(p3) == 5
s3 = sample(Z3, H3, 2^20)
@test length(s3) == 92
@test getRank(s3[[1,1,1,1,1,1,1]], Z3) == 7
ranks = countRanks(s3, Z3)
@test keys(ranks) == Set([7,8,9])
@test [ranks[7], ranks[8], ranks[9]] == [4, 24, 64]

all_zeros = [Bool(0) for k=1:7]
all_ones = [Bool(1) for k=1:7]
P_zero = getPoset(all_zeros, s3, true)
P_one = getPoset(all_ones, s3, areNeighbors)
P_zero_ranks = getPosetRanks(P_zero, s3, Z3)
P_one_ranks = getPosetRanks(P_one, s3, Z3)
printPosetRanks(P_zero_ranks)
printPosetRanks(P_one_ranks)

Z4 = FiniteCyclicGroup([4])
H4 = generateSunletArr(Z4)
s4 = sample(Z4, H4, 2^30)
@test length(s4) == 2328
h4_zeros = [Bool(0) for k=1:13]
h4_ones = [Bool(1) for k=1:13]
h4_P_zero = getPoset(h4_zeros, s4, true)
h4_P_one = getPoset(h4_ones, s4, true)
h4_P_zero_ranks = getPosetRanks(h4_P_zero, s4, Z4)
h4_P_one_ranks = getPosetRanks(h4_P_one, s4, Z4)
printPosetRanks(h4_P_zero_ranks)
printPosetRanks(h4_P_one_ranks)

flips = [[[1],[0],[3]], [[2], [0], [2]], [[3], [0], [1]]]
flipsIndex = [findfirst(x -> x == f, getValidGroupTriples(Z4)) - 1 for f in flips]
flippedRegion = [i in flipsIndex for i=1:13]
flippedMuEta = s4[transpose(flippedRegion)]

containingRegions = []
for k in keys(s4)
    if containsR(k, flippedRegion)
        push!(containingRegions, k)
    end
end
CRsorted = sort!(containingRegions, lt=(x,y)->sum(x)<sum(y))
my_region = CRsorted[1]
my_mueta = s4[my_region]
my_rank = getRank(my_mueta, Z4)

maxrankCR = []
for k in containingRegions
    if getRank(s4[k], Z4) == 16
        push!(maxrankCR, k)
    end
end
MRCRsorted = sort!(maxrankCR, lt=(x,y)->sum(x)<sum(y))
findall(x->x==1, MRCRsorted[4] .- flippedRegion)

altMR = [0 0 1 0 0 0 0 0 1 1 0 0 0]
altME = s4[altMR]

.!([1 1 0] .>= [0 0 1])

function containsR(k, R)
    for i in findall(R)
        if k[i] == 0
            return false
        end
    end
    return true
end