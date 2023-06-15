inlcude("../src/row_reduction.jl")

## get the row reduction for all groups
for group in GROUPS
    row_reduced_csv(group, "data/row_reduced/")
end