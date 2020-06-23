using Test
using PyCall
using BenchmarkTools
using Fermi
using Fermi.Wavefunction
using Fermi.CoupledCluster: RCCSD,RCCD,DFRCCD

psi4.core.be_quiet() #turn off output
psi4.set_num_threads(6)
using LinearAlgebra
BLAS.set_num_threads(6)
# > setup
tol = 1E-14
psi4.set_options(Dict("D_CONVERGENCE" => 14,
                      "E_CONVERGENCE" => 12,
					  "scf_type" => "pk"))
mol2 = psi4.geometry("""
                     O
                     H 1 1.1
                     H 1 1.1 2 104.0
					 symmetry c1
					 """)
e2,wfn2 = psi4.energy("hf/cc-pvtz",mol=mol2,return_wfn=true)
psi4.set_options(Dict("D_CONVERGENCE" => 14,
                      "E_CONVERGENCE" => 12,
                      "scf_type"      => "df"))
e3,wfn3 = psi4.energy("hf/cc-pvtz",mol=mol2,return_wfn=true)
JuWfn3 = Wfn(wfn3)
JuWfn2 = Wfn(wfn2)
printdo=false
#println(@btime RCCD.do_rccd(JuWfn2; doprint=printdo))
println(@btime RCCSD.do_rccsd(JuWfn2; doprint=printdo))
#println(@btime DFRCCD.do_df_rccd(JuWfn3; doprint=printdo))
#println(psi4.energy("ccsd/sto-3g",mol=mol2) - e2)
