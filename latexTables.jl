using LinearAlgebra
include("lambda_cyclic.jl")

function formatTripleLatex(triple)
    return '(' + triple[1][1] + ',' + triple[2][1] + ',' + triple[3][1] + ')'
end

function printTriple(i::Int64, nonTrivialTriples)
    if h == 1
        print("(0,*,*)")
    else
        print(formatTripleLatex(nonTrivialTriples[h-1]))
    end
end

function printHIndex(i::Int64)
    print(" & \$H_{", h, "}\$ & ")
end

function printH(H)
    for j in eachindex(H)
        print(H[j], " & ")
    end
end

function printTableHeader(coords)
    print("\$ (g_1, g_2, g_3) \$ & \$ \\# \$ & ")

    for c in coords
        print("\$", c, "\$ & ")
    end
    
    print("equation \\\\")
end

"""
Given the coordinates and coefficients of a hyperplane, print the corresponding equation
"""
function printEquation(coord, coeff)
    print("\$")
    
    for v in eachindex(coord)
        print(formatVar(coord[v], coeff[v]))
    end

    print(" = 0\$ \\\\")
end

function numToSign(n::Int64)
    if n > 0
        return '+'
    else
        return '-'
    end
end

function formatVar(var, coeff)
    return " \\var " + numToSign(coeff)
end

function getHSupport(H)
    return [i for i in eachindex(H) if H[i] != 0]
end

"""
Print a list of hyperplanes in LaTeX format for copy-and-paste.
"""
function HToLatex(Hs, validTriples, coords)
    nonTrivialTriples = [G for G in validTriples if G[1][1] != 0]
    for h in eachindex(Hs)
        H = Hs[h]
        
        ## (g_1, g_2, g_3)
        printTriple(h, nonTrivialTriples)
        
        ## & $H_{h}$ &
        printHIndex(h)
        
        ## 1 & 0 & ... & 0 &
        printH(H)
        
        ## $\mu_{i} - \eta_{jk} - \eta_{lm} = 0$ \\ 
        H_support = getHSupport(H)
        printEquation(coords[H_support], H[H_support])
        
        println()
    end
end

function generateMuEtaCoordinates(n::Int64)
    mus = [formatMu(i) for i=0:(n-1)]
    etas = [formatEta(i, i-1, n) for i=1:n]
    return vcat(mus, etas)
end

function formatMu(i::Int64)
    return string("mu_{", string(i), "}")
end

function formatEta(i::Int64, j::Int64, n::Int64)
    minIndex = string(min(i % n, j % n))
    maxIndex = string(max(i % n, j % n))
    return string("eta_{", maxIndex, minIndex, "}")
end