using LinearAlgebra

m1 = Dict([4=>[1,2], 5=>[1]])
m2 = Dict([4=>[2], 6=>[1]])
M = [m1, m2]

n = 2

function computeLambda(G, m, l, n)
    func = []
    
    for i in collect(keys(m))
        lower_index = i
        upper_index = (sum([G[j] for j in m[i]]) % n)
        push!(func, (lower_index, upper_index))
    end
    
    return sum([l[I] for I in func])
end

lam = Dict([(5,0)=>-1, (5,1)=>1, (6,0)=>2, (6,1)=>-2, (4,1)=>0, (4,0)=>0])

function getWinner(G, l, n)
    s1 = computeLambda(G, m1, l, n)
    s2 = computeLambda(G, m2, l, n)
    
    if s1 <= s2
        return 1
    else
        return 2
    end
end

## indexes are off, so 
function getVector(g, m, n)
    vector = Array{Int64}(undef, n*5, 1)
    
    g1, g2, g3 = g
    
    a2 = [0 for i=1:n]
    if g2 == n
        a2[1] = 1
    else
        a2[g2 + 1] = 1
    end
    vector[1:n,:] = a2
    
    a3 = [0 for i=1:n]
    if g3 == n
        a3[1] = 1
    else
        a3[g3 + 1] = 1
    end
    vector[n+1:2*n,:] = a3
    
    for i=4:6
        entry = [0 for k=1:n]
        if i in keys(m)
            e = (sum([g[j] for j in m[i]]) % n) + 1
            entry[e] = 1
        end
        vector[(i-2)*n+1:(i-1)*n,:] = entry
    end
    
    return collect(vector)
end

function getMatrix(l, n)
    ## list of triples of group elements
    G = [t for t in collect(Iterators.product(1:n,1:n,1:n)) if (sum(t) % n) == 0]
    
    A = Matrix{Int64}(undef, n*5, length(G))
    
    ## for each group element, compute the winner and then the corresponding vector
    for i=1:length(G)
        g = G[i]
        
        ## computing the winner
        winner = getWinner(g, l, n)
        
        ## getting the corresponding vector
        A[:,i] = getVector(g, M[winner], n)
    end
    
    return A
end

indices(n) = vec([(i,j) for j=0:(n-1), i=4:6])

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

function extended_data(sample, indices, n)
    data = []
    for s in sample
        L = Dict([indices[i]=>s[i] for i=1:3*n])
        B = getMatrix(L, n)
        r = LinearAlgebra.rank(B)
        push!(data, (L, B, r))
    end
    return data
end

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