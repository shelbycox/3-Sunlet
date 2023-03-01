include("groups.jl")
using LinearAlgebra

function generateSunletArr(group::FiniteCyclicGroup)
    validTriples = getValidGroupTriples(group)
    hyperplanes = []
    
    for T in validTriples
        H = getHyperplane(T, group)
        push!(hyperplanes, H)
    end

    return collect(Set(hyperplanes))
end

function getHyperplane(T, group::FiniteCyclicGroup)
    muVector = getMuVector(T, group)
    etaVector = getEtaVector(T, group)

    return [muVector..., etaVector...]
end

function getMuVector(T, group::FiniteCyclicGroup)
    groupElements = getGroupElements(group)
    g1_index = findfirst(x -> x == T[1], groupElements)

    to_return = zeros(Int64, getGroupSize(group))
    to_return[g1_index] = Int(1)
    
    return to_return
end

function getEtaVector(T, group::FiniteCyclicGroup)
    groupElements = getGroupElements(group)
    i = findfirst(x -> x == groupAdd(group, T[1:2]), groupElements)
    j = findfirst(x -> x == T[2], groupElements)

    sign = getEtaSign(i, j)

    return sign*getEtaCoords(i, j, getGroupSize(group))
end

function getEtaSign(i, j)
    if i >= j
        return -1
    else
        return 1
    end
end

function getEtaCoords(i::Int64, j::Int64, n::Int64)
    max_index = max(i,j)
    min_index = min(i,j)

    to_return = zeros(Int64, n-1)
    to_return[min_index : max_index - 1] .= Int(1)

    return to_return
end

##
### This doesn't really need to go here.
##

function getArrangementInequality(v, hyperplanes)
    return [LinearAlgebra.dot(v, H) for H in hyperplanes]
end

function areNeighbors(I, J)
    return sum(I .⊻ J) == 1
end