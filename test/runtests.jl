using Fermi
using Test
using Suppressor

molecules = ["water",   "ammonia",     "benzene", "ethanol", "formaldehyde", "glycine", "methane", "phosphaethene", "acetic_acid", "nitrogen"]
basis =     ["cc-pvtz", "aug-cc-pvdz", "6-31g",   "cc-pvdz", "6-31g*",       "sto-3g",  "cc-pvtz", "3-21g",         "6-311+g",     "cc-pcvdz"]
uhf_boys = [["oh", "cc-pvtz", 0, 2], ["methyl", "cc-pvdz", 0, 2], ["methylene", "cc-pvtz", 0, 3]]
tol = 1E-8

@testset "Fermi" begin
    @time include("test_options.jl")
    @time include("test_output.jl")
    @time include("test_arrays.jl")
    @time include("test_contract.jl")
    @time include("test_molecule.jl")
    @time include("test_diis.jl")
    @time include("test_orbitals.jl")
    @time include("test_phycons.jl")
    @time include("test_integrals.jl")
    @time include("test_RHF.jl")
    @time include("test_UHF.jl")
    @time include("test_MP.jl")
    @time include("test_pT.jl")
    @reset
end
