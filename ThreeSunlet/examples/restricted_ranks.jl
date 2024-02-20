include("../src/row_reduction.jl")

function get_tree_ranks(group)
    L = lambdaGuess(group)
    M = getMatrix(L, group)
    return [LinearAlgebra.rank(M[:, get_tree_i_cols(L, group, i)]) for i=1:2]
end

function print_tree_ranks()
    println("group & \$T_1\$ & \$T_2\$ \\\\")
    println("\\hline")
    for g in GROUPS
        r = get_tree_ranks(g)
        println("\$", string(g), "\$", " & ", r[1], " & ", r[2], " \\\\")
    end
end

function get_g0_ranks(group, g0, k)
    L = lambdaGuess(group)
    M = getMatrix(L, group)
    I = index_triples(group, g0, k)
    return LinearAlgebra.rank(M[:, I])
end

group = FiniteCyclicGroup([5])
get_g0_ranks(group, [2], 1)

function print_triple_g0(g0, k)
    stars = ["*" for i=1:3]
    stars[k] = Base.string(g0)
    return Base.string(stars)
end

function string_g0_ranks_k(group, k)
    to_return = ""
    for g in getGroupElements(group)
        to_return = Base.string(to_return, print_triple_g0(g, k), " & ", get_g0_ranks(group, g, k), " \\\\\n")
    end
    return to_return
end

function string_g0_ranks(group)
    to_return = "triple & rank \\\\\n"
    for k=1:3
        to_return = Base.string(to_return, string_g0_ranks_k(group, k))
    end
    return to_return
end

FILE_PATH = "/Users/shelbycox/Documents/Projects/3-Sunlet/data"

for group in GROUPS
    path = Base.string(FILE_PATH, "/g0_ranks/$(group.structure).txt")
    open(path, "w") do to_write
        write(to_write, string_g0_ranks(group))
    end
end