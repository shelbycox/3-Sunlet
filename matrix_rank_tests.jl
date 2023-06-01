include("lambda_guess.jl")
include("lambda_cyclic.jl")
using CSV, Tables

group = FiniteCyclicGroup([5])
triples = getValidGroupTriples(group)
id = getGroupIdentity(group)

L = lambdaGuess(group)
M = getMatrix(L, group)
r = LinearAlgebra.rank(M)

function my_row_reduction(M::Matrix, group::FiniteCyclicGroup)
    n = getGroupSize(group)
    for i=((n*2)+1):n*3
        M[i,:] = M[i,:] .- M[i-2*n,:]
    end
end

function get_tree_i_cols(L, group::FiniteCyclicGroup, i::Int64)
    triples = getValidGroupTriples(group)
    return [k for k=eachindex(triples) if getLambdaMinimizer(triples[k], L, group) == i]
end

function get_new_col_order(L, group::FiniteCyclicGroup)
    tree1I = get_tree_i_cols(L, group, 1)
    tree2I = get_tree_i_cols(L, group, 2)
    return vcat(tree2I, tree1I)
end

function get_new_row_order(group::FiniteCyclicGroup)
    varOrder = [1, 2, 5, 3, 4]
    n = getGroupSize(group)
    newOrder = []
    for i in varOrder
        newOrder = vcat(newOrder, (n*(i-1)+1):n*i)
    end
    return newOrder
end

function get_row_reduced_M(group::FiniteCyclicGroup)
    L = lambdaGuess(group)
    M = getMatrix(L, group)

    my_row_reduction(M, group)

    newColOrder = get_new_col_order(L, group)
    M = M[:,newColOrder]

    newRowOrder = get_new_row_order(group)
    M = M[newRowOrder,:]

    return M, newColOrder
end

function row_reduced_csv(group::FiniteCyclicGroup, file_path::String)
    M, newColOrder = get_row_reduced_M(group)
    triples = getValidGroupTriples(group)[newColOrder]
    CSV.write("$file_path$(group.structure)_matrix_rr.csv", Tables.table(M), header=[string(T) for T in triples])
end

for group in GROUPS
    row_reduced_csv(group, "data/row_reduced/")
end