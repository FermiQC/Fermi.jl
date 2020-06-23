using Test
using PyCall
#using BenchmarkTools
using Fermi.Wavefunction
#using Fermi.MollerPlesset
using Fermi.Transformation
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
mol2 = psi4.geometry("""
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
e2, wfn2 = psi4.energy("hf/sto-3g", mol = mol2, return_wfn = true)
mints = psi4.core.MintsHelper(wfn2.basisset())
JuWfn2 = Wfn(wfn2)
mol3 = psi4.geometry("""
      1 2
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("reference" => "uhf"))
e3, wfn3 = psi4.energy("hf/sto-3g", mol = mol3, return_wfn = true)
JuWfn3 = Wfn{Float64}(wfn3; unrestricted=true)
@testset "Integral Transformation" begin
    @testset "SmokeRHF" begin
        disk = tei_transform(JuWfn2.uvsr, JuWfn2.Ca, "testdisk")
        mem = tei_transform(JuWfn2.uvsr, JuWfn2.Ca)
        @test disk[:, :, :, :] == mem[:, :, :, :]
    end
    @testset "UHF" begin
        tei_transform(
            JuWfn3.uvsr,
            JuWfn3.Ca,
            JuWfn3.Cb,
            JuWfn3.Ca,
            JuWfn3.Cb,
            "TestTEImixed",
        )
        rhf = tei_transform(
            JuWfn2.uvsr,
            JuWfn2.Ca,
            JuWfn2.Ca,
            JuWfn2.Ca,
            JuWfn2.Ca,
            "TestTEIrhf",
        )
        uhf = tei_transform(
            JuWfn2.uvsr,
            JuWfn2.Ca,
            JuWfn2.Ca,
            JuWfn2.Ca,
            JuWfn2.Ca,
            "TestTEIuhf",
        )
        @test rhf[:, :, :, :] ≈ uhf[:, :, :, :]
    end
    @testset "RHF subset" begin
        _Cao = wfn2.Ca_subset("AO", "OCC")
        _Cav = wfn2.Ca_subset("AO", "VIR")
        subset1 = tei_transform(
            JuWfn2.uvsr,
            JuWfn2.Cav,
            JuWfn2.Cao,
            JuWfn2.Cao,
            JuWfn2.Cav,
            "test",
        )
        subset2 = tei_transform(
            JuWfn2.uvsr,
            JuWfn2.Cav,
            JuWfn2.Cav,
            JuWfn2.Cao,
            JuWfn2.Cao,
            "test2",
        )
        subset3 = tei_transform(
            JuWfn2.uvsr,
            JuWfn2.Cao,
            JuWfn2.Cao,
            JuWfn2.Cao,
            JuWfn2.Cao,
            "test3",
        )
        @test subset1 ≈ mints.mo_eri(_Cav, _Cao, _Cao, _Cav).np
        @test subset2 ≈ mints.mo_eri(_Cav, _Cav, _Cao, _Cao).np
        @test subset3 ≈ mints.mo_eri(_Cao, _Cao, _Cao, _Cao).np
    end
end
