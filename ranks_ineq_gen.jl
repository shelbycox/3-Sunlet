using Test
include("sampling.jl")
include("groups.jl")
include("hyperplanes.jl")
include("lambda_cyclic.jl")

Z3 = FiniteCyclicGroup([3])
H3 = generateSunletArr(Z3)

s3 = sample(Z3, H3, 2^20)
ranks3 = Dict()
for s in keys(s3)
    p = s3[s]
    r = getRank(p, Z3)
    I = getArrangementInequality(p, H3)
    if r in keys(ranks3)
        push!(ranks3[r], I)
    else
        ranks3[r] = [I]
    end
end

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

gen_ranks_ineq(Z3, 2^20)