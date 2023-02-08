using LinearAlgebra, Combinatorics
include("groups.jl")

## choice of basis for the 3-sunlet model
##                      m_1                          m_2
basisExponents = [Dict([4=>[1,2], 5=>[1]]), Dict([4=>[2], 6=>[1]])]

"""
stuff here
groupElementTriple is a list of three group elements.
exponent is the exponent vector (as a list).
lambda is a vector in R^?.
groupStructure is a list of integers that encodes the structure of the (finite) abelian group.
"""
function computeLambda(groupElementTriple, exponent, lambda, group::FiniteCyclicGroup)
    varsToSum = [] ## better name for this!
    
    for i in collect(keys(exponent))
        lower_index = i
        upper_index = groupAdd(group, groupElementTriple[exponent[i]])
        push!(varsToSum, (lower_index, upper_index))
    end
    
    return sum([lambda[I] for I in varsToSum])
end

"""
stuff here
"""
function getLowerLambdaExp(groupElementTriple, lambda, group::FiniteCyclicGroup)
    s1 = computeLambda(groupElementTriple, basisExponents[1], lambda, group)
    s2 = computeLambda(groupElementTriple, basisExponents[2], lambda, group)
    
    if s1 <= s2
        return 1
    else
        return 2
    end
end

"""
stuff here!
## TODO: clean up this code and figure out what it's doing
"""
function getVector(g, basisExponent, group)
    groupSize = getGroupSize(group)
    numFactors = getNumFactors(group)
    vector = Array{Int64}(undef, groupSize*5, 1)
    groupElements = vec(getGroupElements(group))
    
    g1, g2, g3 = g
    
    a2Exponent = [0 for i=1:groupSize] ## initialize a2 exponent indicator vector
    a2Exponent[findfirst(x -> x == g2, groupElements)] = 1
    vector[1:groupSize,:] = a2Exponent
    
    a3Exponent = [0 for i=1:groupSize] ## initialize a3 exponent indicator vector
    a3Exponent[findfirst(x -> x == g3, groupElements)] = 1
    vector[groupSize+1:2*groupSize,:] = a3Exponent
    
    for i=4:6
        entry = [0 for k=1:groupSize] ## initiale a4/5/6 exponent indicator vectors
        if i in keys(basisExponent)
            e = groupAdd(group, [g[j] for j in basisExponent[i]])
            entry[findfirst(x -> x == e, groupElements)] = 1
        end
        vector[(i-2)*groupSize+1:(i-1)*groupSize,:] = entry
    end
    
    return collect(vector)
end

"""
stuff here!
"""
function getMatrix(lambda, group::FiniteCyclicGroup)
    ## list of triples of group elements
    groupElementTriples = getValidGroupTriples(group)
    groupSize = getGroupSize(group)
    A = Matrix{Int64}(undef, groupSize*5, length(groupElementTriples))
    
    ## for each group element, compute the winner and then the corresponding vector
    for i=1:length(groupElementTriples)
        G = groupElementTriples[i]
        
        ## computing the winner
        winner = getLowerLambdaExp(G, lambda, group)
        
        ## getting the corresponding vector
        A[:,i] = getVector(G, basisExponents[winner], group)
    end
    
    return A
end

"""
Given a mu-eta vector, returns a compatible lambda vector.
"""
function to_lambda(mu, eta, group::FiniteCyclicGroup) ## only works for cyclic groups currently! --maybe this same idea will work?
    lambda = Dict()
    groupSize = getGroupSize(group)
    groupElements = getGroupElements(group)
    ## this block should work for any group now
    ## TODO: check that the mus are in the order you think they are wrt group element generation --> they are for Z2 x Z2!
    for i=eachindex(groupElements)
        lambda[(6,groupElements[i])] = maximum([mu[i], 0]) ## set lambda_6^k = mu_k if mu_k > 0, 0 otherwise
        lambda[(5,groupElements[i])] = minimum([mu[i], 0])*(-1) ## set lambda_5^k = -mu_k if mu_k < 0, 0 otherwise
    end

    ## this will now work for every group (assuming coordinates are chosen correctly)
    numFactors = getNumFactors(group)
    lambda[(4,groupElements[1])] = 0 ## set lambda_4^0 = 0
    for i=2:groupSize
        lambda[(4,groupElements[i])] = sum(eta[1:i-1]) ## set lambda_4^k = sum of first k etas
    end

    return lambda
end

function getRank(mu_eta, group::FiniteCyclicGroup)
    groupSize = getGroupSize(group)
    L = to_lambda(mu_eta[1:groupSize], mu_eta[(groupSize+1):end], group)
    return LinearAlgebra.rank(getMatrix(L, group))
end