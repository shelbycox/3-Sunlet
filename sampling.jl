using LinearAlgebra
include("lambda_cyclic.jl")

p(n) = [rand(Float64)-.5 for j=1:n] ## generate a random point

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

function printEquation(coord, coeff)
    print("\$")
    for v in eachindex(coord)
        if v != length(coord)
            if coeff[v+1] > 0
                 print("\\", coord[v], " + ")
            else
                print("\\", coord[v], " - ")
            end
        else
            print("\\", coord[v], " = 0\$ \\\\")
        end
    end
end

function getCoordCoeff(H)
    vars = [i for i in eachindex(H) if H[i] != 0]
    coord = coords[vars]
    coeff = H[vars]
    return coord, coeff
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
        coord, coeff = getCoordCoeff(H)
        printEquation(coord, coeff)
        
        println()
    end
end