include("latexTables.jl")
using Test

@test formatMu(2) == "mu_{2}"
@test formatEta(1,0,3) == "eta_{10}"
@test formatEta(4,2,3) == "eta_{21}"

ME = generateMuEtaCoordinates(3)
@test length(ME) == 6
@test ME == ["mu_{0}", "mu_{1}", "mu_{2}", "eta_{10}", "eta_{21}", "eta_{20}"]

H1 = [1,2,-1,-99999999999999999999,0,0,1,0,2]
H2 = [0,0,0,0,0,0]
H3 = [1,0,0,-1,0,0]
H3_support = getHSupport(H3)
@test getHSupport(H1) == [1,2,3,4,7,9]
@test getHSupport(H2) == []

@test formatVar(ME[1], 1) == "+ \\mu_{0} "
@test formatVar(ME[4], -1) == "- \\eta_{10} "
@test formatVar(ME[2], 0) == "- \\mu_{1} "

@test formatH(H3) == "1 & 0 & 0 & -1 & 0 & 0 & "
@test formatHIndex(10) == " & \$H_{10}\$ & "
