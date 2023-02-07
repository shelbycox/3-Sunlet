using LinearAlgebra, Combinatorics
include("groups.jl")

## choice of basis for the 3-sunlet model
##                      m_1                          m_2
basisExponents = [Dict([4=>[1,2], 5=>[1]]), Dict([4=>[2], 6=>[1]])]

"""
I think I need either a group struct or a setup object.
That way, I can store the group structure, number of factors, group elements, valid group triples, instead of needing to create them every time
"""

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
not sure if I need this anymore
"""
function genData(n, N)
    # for now just recording the rank
    data = []
    for i=1:N
        l = rlam(n)
        A = getMatrix(l, n)
        r = LinearAlgebra.rank(A)
        push!(data, r)
    end
    
    return data
end

## function to generate the indices for Z/nZ
indices(n) = vec([(i,j) for j=0:(n-1), i=4:6])

## helper function for sample(n)
function sample_ranges(n)
    ranges = [0:0]
    for i=1:(n-1)
        push!(ranges, 0:2)
    end

    for i=1:n
        push!(ranges, 0:0)
    end

    for i=1:n
        push!(ranges, -1:1)
    end
    
    return ranges
end

## function to generate the test ranges for lambda
function old_sample(n)
    sample_0 = collect(Iterators.flatten(Iterators.product(sample_ranges(n)...)))
    
    sample = [sample_0[3*n*k+1:3*n*k+3*n] for k=0:-1+Int(trunc(length(sample_0)/(3*n)))]
    
    return sample
end

"""
stuff here!
"""
function get_data(sample, indices, n)
    data = []
    for s in sample
        L = Dict([indices[i]=>s[i] for i=1:3*n])
        B = getMatrix(L, n)
        r = LinearAlgebra.rank(B)
        ## maybe add more to data here?
        push!(data, r)
    end
    return data
end

## compute some statistics of the ranks produced
mma(D) = print(minimum(D), ' ', sum(D)/length(D), ' ', maximum(D));

"""
    ToDo: stuff here
"""
function extended_data(sample, indices, n)
    data = []
    for s in sample
        L = Dict([indices[i]=>s[i] for i=1:3*n])
        B = getMatrix(L, n)
        r = LinearAlgebra.rank(B)
        push!(data, (L, r))
    end
    return data
end

## get the lambda achieving the maximum rank
"""
    ToDo: stuff here
"""
function get_maxes(data, n)
    maxes = []
    for d in data
        if d[3] == 1 + 5*(n-1)
            ## we can just get the matrix again if we want to, so don't record it for now
            push!(maxes, d[1])
        end
    end
    return maxes
end

"""
    ToDo: stuff here
"""
function get_ineq(lambda, n)
    ineq = []
    push!(ineq, ("mu_0 > 0", mu(lambda, 0) > 0))
    for i=1:(n-1)
        for j=1:n-1
            push!(ineq, ("mu_$(i) > eta_$(j)", mu(lambda, i) > eta(lambda, j)))
        end
    end
    return ineq
end

"""
    ToDo: stuff here
"""
function get_gaps(L)
    V = sort(L)
    G = [V[1] - 0.5]
    for i = 1 : (length(V) - 1)
        push!(G, (V[i] + V[i+1])/2)
    end
    push!(G, V[length(V)] + 0.5)
    return G
end

"""
    ToDo: stuff here
"""
function better_sampling(n)
    mus = [] ## store lists of mu_0, mu_1, mu_2, ..., mu_n-1
    ## pick a permutation of the mu_i, assign to 1, ..., n
    for p in collect(Combinatorics.permutations(1:n))
        mu_0 = 2*p[1]

        ## then compute two shifts: mu_0 = -1, = 1
        mu_neg = [2*p[i] - mu_0 - 1 for i=1:n]
        mu_pos = [2*p[i] - mu_0 + 1 for i=1:n]
        
        push!(mus, mu_neg), push!(mus, mu_pos)
    end

    ## mu needs to be paired with eta

    ## now compute some lambda to sample
    test_lam = []
    for mu in mus
        etas = [] ## list to store eta_0, ..., eta_n-1
        gaps = get_gaps(mu)

        eta = [1 for i=1:n]
        push!(etas, [gaps[e] for e in eta])
        while minimum(eta) < n + 1
            ## I need to find the smallest index not maxed out --> j
            j = -1
            for i=1:n
                if eta[i] < n + 1
                    j = i
                    break
                end
            end
            
            ## then reset everything from 1 --> j-1
            for i=1:(j-1)
                eta[i] = 1
            end

            ## and increment the jth index
            eta[j] = eta[j] + 1

            ## then add that eta to the list
            push!(etas, [gaps[e] for e in eta])
        end

        for eta in etas
            lambda = Dict()

            for i=0:(n-1)
                lambda[(6,i)] = mu[i+1] ## set lambda_6^k = mu_k
                lambda[(5,i)] = 0 ## set all lambda_5 = 0
            end

            lambda[(4,0)] = 0 ## set lambda_4^0 = 0
            for i=1:(n-1)
                lambda[(4,i)] = sum(eta[1:i]) ## set lambda_4^k = sum of first k etas
            end

            push!(test_lam, lambda)
        end
    end
    return test_lam
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