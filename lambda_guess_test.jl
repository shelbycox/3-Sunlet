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

m = 100
M = 0
for i=1:100000
    L = lambdaGuess(Z5)
    Mat = getMatrix(L, Z5)
    r = LinearAlgebra.rank(Mat)
    if r < m
        m = r
    elseif r > M
        M = r
    end
end
m, M

L = lambdaGuess(Z5)
M = getMatrix(L, Z5)
r = LinearAlgebra.rank(M)

for i=1:100000
    L = lambdaGuess(Z5)
    newM = getMatrix(L, Z5)
    r = LinearAlgebra.rank(M)
    if M != newM
        print("diff")
        break
    end
end

print(L)

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
    num_rows, num_cols = size(M)
    CSV.write("$file_path$(G.structure)_matrix.csv", Tables.table(M), header=[string(g) for g in getValidGroupTriples(G)])
end

getGuessMatrixCSV(FiniteCyclicGroup([4]), "data/")

G = FiniteCyclicGroup([5])
L = lambdaGuess(G)
M = getMatrix(L, G)
r = LinearAlgebra.rank(M)
num_rows, num_cols = size(M)
target_rank = min(num_rows - 4, num_cols)
CSV.write("Z5_matrix.csv", Tables.table(M), header=[string(g) for g in getValidGroupTriples(G)])

cols = [M[:,i] for i=1:num_cols]
S = collect(subsets(1:num_cols, r))
max_rank_S = [s for s in collect(S) if LinearAlgebra.rank(M[:,s]) == r]


ranks = Set()
for i=1:10000
    L = lambdaGuess(G)
    newM = getMatrix(L, G)
    r = LinearAlgebra.rank(M)
    push!(ranks, r)
    if r != 15
        print(i)
        break
    end
end
ranks