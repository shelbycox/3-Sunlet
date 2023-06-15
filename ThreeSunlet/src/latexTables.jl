using LinearAlgebra
include("lambda_cyclic.jl")

"""
Print a list of hyperplanes in LaTeX format for copy-and-paste.
"""
function HToLatex(Hs, group::FiniteCyclicGroup)
    nonTrivialTriples = [(0,0,0), [G for G in getValidGroupTriples(group) if G[1] != getGroupIdentity(group)]...]
    coords = generateMuEtaCoordinates(group)

    to_return = formatTableHeader(coords)
    
    for h in eachindex(Hs) ## TODO: indices of Hs and nonTrivialTriples are not the same!
        H = Hs[h]
        
        ## (g_1, g_2, g_3)
        to_return = string(to_return, formatTriple(nonTrivialTriples[h], group))

        ## & $H_{h}$ &
        to_return = string(to_return, formatHIndex(h))
        
        ## 1 & 0 & ... & 0 &
        to_return = string(to_return, formatH(H))
        
        ## $\mu_{i} - \eta_{jk} - \eta_{lm} = 0$ \\ 
        H_support = getHSupport(H)
        to_return = string(to_return, formatEquation(coords[H_support], H[H_support]), "\n")
    end
    
    to_return = string(to_return, "\\hline \n")

    return to_return
end

function formatTableHeader(coords)
    to_return = "\$ (g_1, g_2, g_3) \$ & \$ \\# \$ & "

    for c in coords
        to_return = string(to_return, "\$ \\", c, "\$ & ")
    end
    
    return string(to_return, "equation \\\\ \n \\hline \n")
end

function formatTriple(triple, group::FiniteCyclicGroup)
    if triple[1] == getGroupIdentity(group)
        return "(0,*,*)"
    else
        return string('(', triple[1], ',', triple[2], ',', triple[3], ')')
    end
end

function formatHIndex(i::Int64)
    return string(" & \$H_{", i, "}\$ & ")
end

function formatH(H)
    return string([string(H[j], " & ") for j in eachindex(H)]...)
end

"""
Given the coordinates and coefficients of a hyperplane, print the corresponding equation
"""
function formatEquation(coord, coeff)
    to_return = "\$"
    
    for v in eachindex(coord)
        to_return = string(to_return, formatVar(coord[v], coeff[v]))
    end

    return string(to_return, "= 0\$ \\\\")
end

function numToSign(n::Int64)
    if n > 0
        return '+'
    else
        return '-'
    end
end

function formatVar(var, coeff)
    return string(numToSign(coeff), " \\", var, " ")
end

function getHSupport(H)
    return [i for i in eachindex(H) if H[i] != 0]
end

function generateMuEtaCoordinates(n::Int64)
    mus = [formatMu(i) for i=0:(n-1)]
    etas = [formatEta(i, i-1, n) for i=1:n-1]
    return vcat(mus, etas)
end

function generateMuEtaCoordinates(group::FiniteCyclicGroup)
    return generateMuEtaCoordinates(getGroupSize(group))
end

function formatMu(i::Int64)
    return string("mu_{", string(i), "}")
end

function formatEta(i::Int64, j::Int64, n::Int64)
    minIndex = string(min(i % n, j % n))
    maxIndex = string(max(i % n, j % n))
    return string("eta_{", maxIndex, minIndex, "}")
end