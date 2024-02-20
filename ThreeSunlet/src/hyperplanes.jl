include("groups.jl")
using LinearAlgebra

function generateSunletArr(group::FiniteCyclicGroup)
    validTriples = getValidGroupTriples(group)
    zeroHyperplane = [Int(i==1) for i=1:2*getGroupSize(group)-1]
    hyperplanes = [zeroHyperplane]
    
    for T in validTriples
        if T[1] != getGroupIdentity(group) 
            H = getHyperplane(T, group)
            push!(hyperplanes, H)
        end
    end

    return hyperplanes
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

## FIXME: make this work for non-cyclic groups!
function getEtaCoords(i::Int64, j::Int64, n::Int64)
    max_index = max(i,j)
    min_index = min(i,j)

    to_return = zeros(Int64, n-1)
    to_return[min_index : max_index - 1] .= Int(1)

    return to_return
end


"""
Works for mu-eta coordinates.
"""
function getArrangementInequality(v, hyperplanes)
    return [LinearAlgebra.dot(v, H) > 0 for H in hyperplanes]
end

function areNeighbors(I, J)
    return sum(I .⊻ J) == 1
end

function areM0Neighbors(I, J)
    return sum(I[2:end] .⊻ J[2:end]) == 1
end

function containsNeighbor(I, Js)
    for J in Js
        if areNeighbors(I, J)
            return true
        end
    end
    return false
end

## TODO: test this
function distToRegion(source, target, regions)
    count = 0
    curr = [source]
    while target ∉ curr
        count = count + 1
        new = [R for R in regions if containsNeighbor(R, curr)]
        curr = new
    end
    return count
end

function lambda_to_mueta(lambda, group::FiniteCyclicGroup)
    groupSize = getGroupSize(group)
    elements = getGroupElements(group)
    eta = [lambda[(4, elements[i+1])] - lambda[(4, elements[i])] for i=1:groupSize-1]
    mu = [lambda[(6, elements[i])] - lambda[(5, elements[i])] for i=1:groupSize]
    return vcat(mu, eta)
end