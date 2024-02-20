include("sampling.jl")
include("groups.jl")
include("hyperplanes.jl")
include("lambda_cyclic.jl")

function gen_ranks_ineq(group::FiniteCyclicGroup, N::Int64)
    H = generateSunletArr(group)
    S = sample(group, H, N)
    ranks = Dict()
    for s in keys(S)
        p = S[s]
        r = getRank(p, group)
        I = getArrangementInequality(p, H)
        if r in keys(ranks)
            push!(ranks[r], I)
        else
            ranks[r] = [I]
        end
    end

    return ranks
end