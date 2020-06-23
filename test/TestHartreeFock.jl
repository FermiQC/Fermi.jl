using Test
using Fermi.HartreeFock
using Fermi
#psi4 = pyimport("psi4")
#psi4.core.be_quiet() #turn off output
mol = psi4.geometry("""
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type" => "pk", "d_convergence" => 14))
#e, wfn = psi4.energy("hf/sto-3g", mol = mol2, return_wfn = true)
wfn = RHFWfn(mol)
println("Created wfn")
RHFCompute(wfn,doprint=true)
#@testset "CoupledCluster" begin
#    @testset "Smoke" begin
#        @test do_rccd(JuWfn2, 40, false) ≈ -0.07015050066089029
#        @test do_rccd(JuWfn3, 40, false) ≈ -0.07015050066089029
#    end
#end
