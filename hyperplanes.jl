include("groups.jl")

function generateSunletArr(group::FiniteCyclicGroup)
    groupSize = prod(group.structure)
    validTriples = getValidGroupTriples(group)
    groupElements = vec(getGroupElements(group))
    numFactors = getNumFactors(group)
    
    hyperplanes = []
    
    for T in validTriples
        muVector = getMuVector(T, group)
        etaVector = getEtaVector(T, group)

        H = [muVector..., etaVector...]
        
        push!(hyperplanes, H)
    end
    return List(Set(hyperplanes))
end

function getMuVector(T, group::FiniteCyclicGroup)
    to_return = zeros(getGroupSize(group))
    to_return[T[1]] = 1
    
    return to_return
end

function getEtaVector(T, group::FiniteCyclicGroup)
    i = groupAdd(group, T[1:2])
    j = T[2]

    sign = getEtaSign(i, j)

    return sign.*getEtaCoords(i, j, getGroupSize(group))
end

function getEtaSign(i, j)
    if i >= j
        return 1
    else
        return -1
    end
end

function getEtaCoords(i::Int64, j::Int64, n::Int64)
    max_index = max(i,j)
    min_index = min(i,j)

    to_return = zeros(n)
    to_return[min_index:max_index] = 1

    return to_return
end