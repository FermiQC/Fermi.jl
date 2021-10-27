using Fermi
using Test
using Suppressor

molecules = ["water",   "ammonia",     "benzene", "ethanol", "formaldehyde", "glycine", "methane", "phosphaethene", "acetic_acid", "nitrogen"]
basis =     ["cc-pvtz", "aug-cc-pvdz", "6-31g",   "cc-pvdz", "6-31g*",       "sto-3g",  "cc-pvtz", "3-21g",         "6-311+g",     "cc-pcvdz"]
uhf_boys = [["oh", "cc-pvtz", 0, 2], ["methyl", "cc-pvdz", 0, 2], ["methylene", "cc-pvtz", 0, 3]]
tol = 1E-8

@testset "Fermi" begin
    include("test_options.jl")
    include("test_output.jl")
    include("test_arrays.jl")
    include("test_contract.jl")
    include("test_geometry.jl")
    include("test_diis.jl")
    include("test_orbitals.jl")
    include("test_phycons.jl")
    include("test_RHF.jl")
    include("test_UHF.jl")
    include("test_MP.jl")
    include("test_pT.jl")
    @reset
end
