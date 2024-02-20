include("groups.jl")

function containsInv(group::FiniteCyclicGroup, X, g)
    e = getGroupIdentity(group)
    for x in X
        if groupAdd(group, [x, g]) == e
            return true
        end
    end
    return false
end

function getNonInvSet(group::FiniteCyclicGroup)
    X = []
    ## loop through the group elements
    elements = Vector(getGroupElements(group))
    e = getGroupIdentity(group)
    for i=eachindex(elements)
        g = elements[i]
        if !(g == e) & !(containsInv(group, X, g))
            push!(X, g)
        end
    end
    return X
end

function lambdaGuess(group::FiniteCyclicGroup)
    X = getNonInvSet(group)
    groupElements = getGroupElements(group)
    e = getGroupIdentity(group)
    lambda = Dict()
    mu = Dict()

    ## pick lambda_g random, small, and positive
    for g in groupElements
        lambda[g] = rand()
    end
    lambda[e] = -10

    ## for each g \in X, pick mu_g so that 
    ##lambda_e - lambda_-g < mu_g < lambda_(g+h) - h
    ##for all h in G
    ## mu_g should end up small and negative... I think
    for g in X
        max = minimum([lambda[groupAdd(group, [g, h])] - lambda[h] for h in groupElements])
        min = lambda[e] - lambda[getInverse(group, g)]
        mu[g] = -5
    end
    mu[e] = 1

    ## set mu_-g = -mu_g
    finishMuGuess(group, mu, X)

    return compileLambda(lambda, mu)
end

function compileLambda(L, M)
    lambda = Dict()
    for k in keys(L)
        lambda[(5, k)] = 0
        lambda[(6, k)] = M[k]
        lambda[(4, k)] = L[k]
    end

    return lambda
end

function finishMuGuess(group::FiniteCyclicGroup, partialMu, X)
    elements = getGroupElements(group)
    for g in elements
        if g ∉ X
            partialMu[g] = -partialMu[getInverse(group, g)]
        end
    end
    return partialMu
end