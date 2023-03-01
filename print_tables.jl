include("groups.jl")
include("hyperplanes.jl")
include("latexTables.jl")

Z3 = FiniteCyclicGroup([3])
H3 = generateSunletArr(Z3)
open("H3_table.txt", "w") do to_write
    write(to_write, HToLatex(H3, Z3))
end

Z4 = FiniteCyclicGroup([4])
H4 = generateSunletArr(Z4)
open("H4_table.txt", "w") do to_write
    write(to_write, HToLatex(H4, Z4))
end

Z2Z2 = FiniteCyclicGroup([2,2])
H22 = generateSunletArr(Z2Z2)
open("H22_table.txt", "w") do to_write
    write(to_write, HToLatex(H22, Z2Z2))
end