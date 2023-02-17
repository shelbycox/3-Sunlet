include("latexTables.jl")
using Test

@test formatMu(2) == "mu_{2}"
@test formatEta(1,0,3) == "eta_{10}"
@test formatEta(4,2,3) == "eta_{21}"

ME = generateMuEtaCoordinates(3)
@test length(ME) == 6
@test ME == ["mu_{0}", "mu_{1}", "mu_{2}", "eta_{10}", "eta_{21}", "eta_{20}"]

