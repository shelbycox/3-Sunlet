include("cyclic.jl")

indices_6 = indices(6)
sample_6 = collect(Iterators.product(0, 0:2, 0:2, 0:2, 0:2, 0:2, 0, 0, 0, 0, 0, 0, -1:1, -1:1, -1:1, -1:1, -1:1, -1:1));

D6 = extended_data(sample_6, indices_6, 6)
M6 = [d[1] for d in D6 if d[3] == 26]

indices_5 = indices(5)
sample_5 = collect(Iterators.product(0, 0:2, 0:2, 0:2, 0:2, 0, 0, 0, 0, 0, -1:1, -1:1, -1:1, -1:1, -1:1));

D5 = extended_data(sample_5, indices_5, 5)
M5 = [d[1] for d in D5 if d[3] == 21]

indices_4 = indices(4)
sample_4 = collect(Iterators.product(0, 0:2, 0:2, 0:2, 0, 0, 0, 0, -1:1, -1:1, -1:1, -1:1));

D4 = extended_data(sample_4, indices_4, 4)
M4 = [d[1] for d in D4 if d[3] == 16]