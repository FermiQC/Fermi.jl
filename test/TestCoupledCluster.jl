using Test
using PyCall
using Fermi
using Fermi.Wavefunction
using Fermi.CoupledCluster: RCCSD, RCCD, ROCCD, DFRCCD, UCCSD, mRCCD, mRCCSD


psi4.core.be_quiet() #turn off output
psi4.set_num_threads(6)
tol = 1E-14
mol2 = psi4.geometry("""
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type" => "pk", "d_convergence" => 14))
e, wfn2 = psi4.energy("hf/sto-3g", mol = mol2, return_wfn = true)
psi4.set_options(Dict("scf_type"=>"df"))
e, wfn2_df = psi4.energy("hf/sto-3g", mol = mol2, return_wfn = true)
psi4.set_options(Dict("scf_type"=>"pk"))
JuWfn2 = Wfn(wfn2)
JuWfn2_s = Wfn{Float32}(wfn2)
JuWfn2_df = Wfn(wfn2_df)
mol3 = psi4.geometry("""
      1 2
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("reference" => "rohf"))
e, wfn3 = psi4.energy("hf/sto-3g",mol=mol3,return_wfn=true)
JuWfn3 = Wfn(wfn3)
psi4.set_options(Dict("reference" => "uhf", "scf_type" => "pk"))
e,wfn4 = psi4.energy("hf/sto-3g",mol=mol3,return_wfn=true)
JuWfn4 = Wfn{Float64}(wfn4; unrestricted=true)
#JuWfn3 = Wfn(wfn2, Float64, true, true) #disk based CCD is currently NOT working
@testset "CoupledCluster" begin
    @testset "Smoke" begin
        #@test RCCD.do_rccd(JuWfn2) ≈ -0.07015050066089029
        #@test mRCCD.do_rccd(JuWfn2) ≈ -0.07015050066089029
        #@test RCCD.do_rccd(JuWfn2_s) ≈ -0.0701504929121782
        #@test RCCSD.do_rccsd(JuWfn2) ≈ -0.070680102078571
        @test mRCCSD.do_rccsd(JuWfn2) ≈ -0.070680102078571
        #@test RCCSD.do_rccsd(JuWfn2_s) ≈ -0.07068009398829002
        #E = UCCSD.do_uccsd(JuWfn2; doprint=true)
        #println(DFRCCD.do_df_rccd(JuWfn2_df; doprint=true))
        #ROCCD.do_roccd(JuWfn3, 40, doprint=true)
    end
end
