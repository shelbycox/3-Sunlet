using LinearAlgebra

p(n) = [rand(Float64)-.5 for j=1:n] ## get a random point

"""
Returns a list of hyperplanes in the arrangement.
"""
function getH(k, full=true)
    H = []
    if full
        for i=1:(2*k)
            v = zero(rand(2*k))
            v[i] = 1
            push!(H, v)
        end
        
        for i=1:(2*k)
            for j=(i+1):(2*k)
                v = zero(rand(2*k))
                v[i] = 1
                v[j] = -1
                push!(H, v)
            end
        end
    else
        v = zero(rand(2*k))
        v[1] = 1
        push!(H, v)
        for i=1:k
            for j=1:k
                v = zero(rand(2*k))
                v[i] = 1
                v[k + j] = -1
                push!(H, v)
            end
        end
    end
    return H
end

"""
Reflects a vector $v$ over the hyperplane with normal vector $a$
"""
function reflect(v, a)
    return v .- a*2*(LinearAlgebra.dot(v,a))/LinearAlgebra.dot(a,a)
end

"""
Returns the orbit of v under the hyperplane arrangment H.
"""
function getOrbit(v, H)
    k = Int(length(v)/2)
    orbit = [v]
    ineqs = [ineq(v, k)]
    to_check = [v]
    while length(to_check) > 0
        new_to_check = []
        for u in to_check
            for h in H
                r = reflect(u, h)
                I = ineq(r, k)
                if I ∉ ineqs
                    push!(orbit, r)
                    push!(new_to_check, r)
                    push!(ineqs, I)
                end
            end
        end
        to_check = new_to_check
    end
    return orbit
end

"""
Returns a matrix that records the Weyl chamber of v.
"""
function ineq(v, k, full=false)
    if full
        return [[v[i] > 0 for i=1:(2*k)] [v[i] > v[j] for i=1:(2*k), j=1:(2*k)]]
    end
    top_row = [1 < 0 for i=1:k]
    top_row[1] = (v[1] > 0)
    return [top_row [v[i] > v[j] for i=1:k, j=(k+1):(2*k)]]
end

"""
Project the vector v onto the hyperplane with normal vector a.
"""
function project(v, a)
    return v - a*(LinearAlgebra.dot(v, a))/(LinearAlgebra.dot(a,a))
end

"""
Project a set of vectors orbit onto the hyperplane sum(eta_j) = 0.
"""
function getProj(orbit, k)
    z = zero(rand(2*k))
    for i=(k+1):(2*k)
        z[i] = 1
    end
    
    P = [project(u, z) for u in orbit]
    
    sample = []
    Is = []
    for u in P
        I = ineq(u, k)
        if I ∉ Is
            push!(sample, u)
            push!(Is, I)
        end
    end
    
    return sample
end

"""
Checks if u, v are in neighboring Weyl chambers.
"""
function neighbors(orbit, v, k, full=false)
    return [u for u in orbit if isNeighbor(u, v, k, full)]
end

"""
Returns a list of all the neighbors of v in orbit.
"""
function isNeighbor(u, v, k, full=false)
    return sum(ineq(u, k, full) .⊻ ineq(v, k, full)) == 1
end