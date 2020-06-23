using Test
using PyCall
#using BenchmarkTools
using Fermi.Wavefunction
using Fermi.MollerPlesset
using Fermi.CoupledCluster: UCCSD
using Fermi
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
psi4.set_options(Dict("D_CONVERGENCE" => 10, "scf_type" => "pk"))
mol2 = psi4.geometry("""
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
e2, wfn2 = psi4.energy("hf/sto-3g", mol = mol2, return_wfn = true)
psi4.set_options(Dict("scf_type" => "df"))
e4,wfn4 = psi4.energy("hf/sto-3g",mol=mol2,return_wfn=true)
JuWfn4 = Wfn(wfn4)
psi4.set_options(Dict("scf_type" => "pk"))
mol3 = psi4.geometry("""
      1 2
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("reference" => "uhf"))
e3, wfn3 = psi4.energy("hf/sto-3g", mol = mol3, return_wfn = true)
JuWfn2 = Wfn{Float64}(wfn2; unrestricted=true)
JuWfn3 = Wfn{Float64}(wfn3; unrestricted=true)
@testset "MP2" begin
    @test do_rmp2(JuWfn2) ≈ -0.04914964480386458
    @test do_ump2(JuWfn3) ≈ -0.03588729625230
    @test do_rmp2(JuWfn2) ≈ do_ump2(JuWfn2)
    @test do_direct_rmp2(JuWfn2) ≈ -0.04914964480386458
    @test do_df_rmp2(JuWfn4) ≈ -0.04913505451294127
end
