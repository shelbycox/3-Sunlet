struct FiniteCyclicGroup
    structure::Vector{Int64}
end

function getGroupSize(group::FiniteCyclicGroup)
    return prod(group.structure)
end

function getNumFactors(group::FiniteCyclicGroup)
    return length(group.structure)
end

function getGroupIdentity(group)
    return zeros(Int64, getNumFactors(group))
end

function groupAdd(group::FiniteCyclicGroup, elementsToAdd)
    return sum(elementsToAdd) .% group.structure
end

function isSumZeroTriple(threeElements, group::FiniteCyclicGroup)
    return groupAdd(group, threeElements) == getGroupIdentity(group)
end

"""
Given the structure of a finite abelian group, returns a list of all elements (each element is encoded as a list).
"""
function getGroupElements(group::FiniteCyclicGroup)
    groupElementForm = [0:j-1 for j in group.structure]
    return vec([collect(g) for g in collect(Iterators.product(groupElementForm...))])
end

"""
Given the structure of a finite abelian group, returns a list of all triples of group elements that sum to zero.
"""
function getValidGroupTriples(group::FiniteCyclicGroup)
    groupElements = getGroupElements(group)
    groupElementTriples = [collect(T) for T in collect(Iterators.product(groupElements, groupElements, groupElements))]
    sumZeroTriples = [T for T in groupElementTriples if groupAdd(group, T) == zeros(getNumFactors(group))]
    return sumZeroTriples
end