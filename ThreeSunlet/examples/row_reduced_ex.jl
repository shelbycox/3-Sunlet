include("../src/row_reduction.jl")

## get the row reduction for all groups
for group in GROUPS
    row_reduced_csv(group, "data/row_reduced/")
end

for group in GROUPS
    L = lambdaGuess(group)
    M = get_row_reduced_M(group)[1]
    n = getGroupSize(group)
    T2_end = length(get_tree_i_cols(L, group, 2))
    println(group.structure, " ", LinearAlgebra.rank(M[1:3*n, 1:T2_end]), " ", LinearAlgebra.rank(M[(3*n+1):end,(T2_end + 1):end]))
end

group = FiniteCyclicGroup([6])
M = get_row_reduced_M(group)[1]
M[:, 1]