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

"""
Given a group element g, find its (additive) inverse in the group
"""
function getInverse(group::FiniteCyclicGroup, g)
    return [(-g[i][1] + group.structure[i]) % group.structure[i] for i=eachindex(group.structure)]
end

GROUP_STRUCTURES = [[3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14], [15], [16], [17], [18], [2,2], [2,2,2], [2,2,2,2], [2,2,2,2,2], [4,2,2], [3,3], [3,3,3]]
GROUPS = [FiniteCyclicGroup(s) for s in GROUP_STRUCTURES]