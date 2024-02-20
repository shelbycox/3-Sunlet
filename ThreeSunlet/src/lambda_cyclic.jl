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
    varsToSum = collectVarsToSum(group, exponent, groupElementTriple)
    return sum([lambda[I] for I in varsToSum])
end

function collectVarsToSum(group::FiniteCyclicGroup, exponent, groupElementTriple)
    to_return = [] ## better name for this!
    
    for i in collect(keys(exponent))
        upper_index = groupAdd(group, groupElementTriple[exponent[i]])
        push!(to_return, (i, upper_index))
    end

    return to_return
end

"""
stuff here
"""
function getLambdaMinimizer(groupElementTriple, lambda, group::FiniteCyclicGroup)
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
basis element is m1 or m2
g is a triple of group elements
group is the finite abelian group
"""
function getVector(g, basisExponent, group)
    groupSize = getGroupSize(group)
    numFactors = getNumFactors(group)
    vector = Array{Int64}(undef, groupSize*6, 1)
    groupElements = vec(getGroupElements(group))
    
    g1, g2, g3 = g

    a1Exponent = [0 for i=1:groupSize]
    a1Exponent[findfirst(x -> x == g1, groupElements)] = 1
    vector[1:groupSize,:] = a1Exponent
    
    a2Exponent = [0 for i=1:groupSize] ## initialize a2 exponent indicator vector
    a2Exponent[findfirst(x -> x == g2, groupElements)] = 1
    vector[groupSize+1:2*groupSize,:] = a2Exponent
    
    a3Exponent = [0 for i=1:groupSize] ## initialize a3 exponent indicator vector
    a3Exponent[findfirst(x -> x == g3, groupElements)] = 1
    vector[2*groupSize+1:3*groupSize,:] = a3Exponent
    
    for i=4:6
        entry = [0 for k=1:groupSize] ## initiale a4/5/6 exponent indicator vectors
        if i in keys(basisExponent)
            e = groupAdd(group, [g[j] for j in basisExponent[i]])
            entry[findfirst(x -> x == e, groupElements)] = 1
        end
        vector[(i-1)*groupSize+1:i*groupSize,:] = entry
    end
    
    return collect(vector)
end

"""
stuff here!
"""
function getMatrix(lambda, group::FiniteCyclicGroup)
    groupElementTriples = getValidGroupTriples(group)
    groupSize = getGroupSize(group)
    A = Matrix{Int64}(undef, groupSize*6, length(groupElementTriples))
    
    for i=eachindex(groupElementTriples)
        G = groupElementTriples[i]
        minimizer = getLambdaMinimizer(G, lambda, group)
        A[:,i] = getVector(G, basisExponents[minimizer], group)
    end
    
    return A
end

"""
Given a mu-eta vector, returns a compatible lambda vector.
"""
function to_lambda(mu, eta, group::FiniteCyclicGroup) ## only works for cyclic groups currently! --maybe this same idea will work?
    groupElements = getGroupElements(group)
    
    muLambda = muToLambda(mu, groupElements)
    etaLambda = etaToLambda(eta, groupElements)
    
    return merge(muLambda, etaLambda)
end

function muToLambda(mu, groupElements)
    to_return = Dict()

    for i=eachindex(groupElements)
        to_return[(6,groupElements[i])] = maximum([mu[i], 0]) ## set lambda_6^k = mu_k if mu_k > 0, 0 otherwise
        to_return[(5,groupElements[i])] = minimum([mu[i], 0])*(-1) ## set lambda_5^k = -mu_k if mu_k < 0, 0 otherwise
    end

    return to_return
end

function etaToLambda(eta, groupElements)
    to_return = Dict()

    to_return[(4,groupElements[1])] = 0 ## set lambda_4^0 = 0
    for i=eachindex(groupElements[2:end])
        ## TODO: check for bugs here
        to_return[(4,groupElements[i+1])] = sum(eta[1:i]) ## set lambda_4^k = sum of first k etas
    end

    return to_return
end

function getRank(mu_eta, group::FiniteCyclicGroup)
    groupSize = getGroupSize(group)
    L = to_lambda(mu_eta[1 : groupSize], mu_eta[(groupSize + 1):end], group)
    return LinearAlgebra.rank(getMatrix(L, group))
end