using Test
using IterTools
using CSV, Tables
include("lambda_guess.jl")
include("lambda_cyclic.jl")

Z5 = FiniteCyclicGroup([5])
E5 = getGroupElements(Z5)
g1 = [1]
g2 = [2]
g3 = [4]
X = [g2]
Y = [g3]
Z = [g2, g3]

@test containsInv(Z5, X, g1) == false
@test containsInv(Z5, Y, g1) == true
@test containsInv(Z5, Z, g1) == true

N = getNonInvSet(Z5)
@test length(N) == 2

L = lambdaGuess(Z5)
M = getMatrix(L, Z5)
r = LinearAlgebra.rank(M)

function latexMatrix(matrix, file_name, groupTriples)
    # open(file_name, "w")
    
    print(file_name, "\\begin{array}{*{", size(matrix)[2], "}c}", "\n")

    for i=eachindex(groupTriples)
        t = groupTriples[i]
        print(file_name, "\\begin{turn}{90} ", t, "\\end{turn}")
        if i == length(groupTriples)
            print(file_name, " \\\\ \n")
        else
            print(file_name, " & ")
        end
    end

    for i=1:size(matrix)[1]
        for j=1:size(matrix)[2]
            print(file_name, matrix[i,j])

            if j != size(matrix)[2]
                print(file_name, " & ")
            end
        end
        if i != size(matrix)[2]
            print(file_name, " \\\\ \n")
        else
            print(file_name, "\n")
        end
    end

    print(file_name, "\\end{array}")
end

gt = getValidGroupTriples(Z5)
latexMatrix(M, stdout, gt)

for k in keys(L)
    println("\\lambda_", k[1], "^", k[2][1], " &= ", L[k], " \\\\")
end

function getGuessMatrixCSV(group::FiniteCyclicGroup, file_path::String)
    G = group
    L = lambdaGuess(G)
    M = getMatrix(L, G)
    r = LinearAlgebra.rank(M)
    CSV.write("$file_path$(G.structure)_matrix.csv", Tables.table(M), header=[string(g) for g in getValidGroupTriples(G)])
end

getGuessMatrixCSV(FiniteCyclicGroup([4]), "data/")