struct FiniteCyclicGroup
    structure::Vector{Int64}
    # groupSize::Int64
    # numFactors::Int64
    # elements::Any
    # validElementTriples::Any
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
Given the structure of a finite abelian group, returns a list of all elements (each element is encoded as a list).
"""
function getGroupElements(group::FiniteCyclicGroup)
    groupElementForm = [0:j-1 for j in group.structure]
    return [collect(g) for g in collect(Iterators.product(groupElementForm...))]
end

"""
Computes the size of a finite cyclic group.
"""
function getGroupSize(group::FiniteCyclicGroup)
    return prod(group.structure)
end

function getNumFactors(group::FiniteCyclicGroup)
    return length(group.structure)
end

"""
Sums a list of group elements.
"""
function groupAdd(group::FiniteCyclicGroup, elementsToAdd)
    return sum(elementsToAdd) .% group.structure
end