using Fermi
using Test
using Suppressor

molecules = ["water",   "ammonia",     "benzene", "ethanol", "formaldehyde", "glycine", "methane", "phosphaethene", "acetic_acid", "nitrogen"]
basis =     ["cc-pvtz", "aug-cc-pvdz", "6-31g",   "cc-pvdz", "6-31g*",       "sto-3g",  "cc-pvtz", "3-21g",         "6-311+g",     "cc-pcvdz"]

tol = 1E-8
@testset "Fermi" begin
    include("test_options.jl")
    #include("test_output.jl")
    #include("test_RHF.jl")
    include("test_MP.jl")
    #include("test_CCSD.jl")
    include("test_pT.jl")
    @reset
end
