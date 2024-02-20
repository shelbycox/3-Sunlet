include("../src/groups.jl")
include("../src/lambda_cyclic.jl")
include("../src/sampling.jl")

N = 2^28
max_ranks = Dict()

for g in GROUPS
    ## sample until finding something of maximal rank or hitting 2^28 samples
    max_ranks[g] = 0
    n = getGroupSize(g)
    for i=1:N
        r = getRank(getSamplePoint(g), g)
        if r > max_ranks[g]
            max_ranks[g] = r
            if r == min(n^2, 5*(n-1) + 1)
                break
            end
        end
    end
end

max_ranks