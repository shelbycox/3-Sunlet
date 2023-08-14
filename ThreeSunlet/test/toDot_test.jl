include("../src/toDot.jl")

group = FiniteCyclicGroup([3])
P = gen_graph_poset(group, 2^10)