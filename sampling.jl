include("hyperplanes.jl")

function sample(group::FiniteCyclicGroup, H, N::Int)
    sample = Dict()
    for i=1:N
        new_point = getSamplePoint(group)
        sample[getArrangementInequality(new_point, H)] = new_point
    end
    return sample
end

function getSamplePoint(group::FiniteCyclicGroup)
    return [rand(Float64) - 0.5 for k=1:2*getGroupSize(group)-1]
end

function countRanks(sample, group::FiniteCyclicGroup)
    ranks = [getRank(s, group) for s in collect(values(sample))]
    R = sort!(collect(Set(ranks)))
    return Dict([k=>sum([r == k for r in ranks]) for k in R])
end

function getNeighbors(I, sample, neighbors)
    return [J for J in collect(keys(sample)) if neighbors(I, J)]    
end

function getPoset(baseRegion, sample, toggleM0)
    poset = Dict(1=>[baseRegion]) ## start with the parochial chamber as a base region
    visited = [baseRegion]
    check = 0
    while length(poset) > check
        check = length(poset)
        next_level = []
        for I in poset[check]
            if toggleM0
                neighbors = getNeighbors(I, sample, areM0Neighbors)
            else
                neighbors = getNeighbors(I, sample, areNeighbors)
            end
            
            for N in neighbors
                if N ∉ visited
                    push!(next_level, N)
                    push!(visited, N)
                end
            end
        end
        if length(next_level) > 0
            poset[check + 1] = next_level
        end
    end
    return poset
end

function getPosetRanks(poset, sample, group::FiniteCyclicGroup)
    return Dict([k=>Set([getRank(sample[I], group) for I in poset[k]]) for k in collect(keys(poset))])
end

function printPosetRanks(poset_ranks)
    for level in sort!(collect(keys(poset_ranks)))
        println(level, " ", poset_ranks[level])
    end
end