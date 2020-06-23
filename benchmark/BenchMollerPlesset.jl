using Test
using PyCall
using BenchmarkTools
using Fermi.Wavefunction
using Fermi.MollerPlesset
using LinearAlgebra.BLAS
using Fermi
BLAS.set_num_threads(4)

psi4.core.set_num_threads(4)
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
psi4.set_options(Dict("D_CONVERGENCE" => 10,
					  "scf_type" => "pk"))
mol2 = psi4.geometry("""
                     pubchem:ethane
					 symmetry c1
					 """)
e2,wfn2 = psi4.energy("hf/cc-pvtz",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2,Float64,false,false)
#mol3 = psi4.geometry("""
#					 1 2
#					 O
#					 H 1 1.1
#					 H 1 1.1 2 104.0
#					 symmetry c1
#					 """)
#psi4.set_options(Dict("reference" => "uhf"))
#e3,wfn3 = psi4.energy("hf/cc-pvdz",mol=mol3,return_wfn=true)
println("Starting benchmark!")
#println(@benchmark Wfn(wfn2,Float64,false,false))
println(@benchmark Wfn(wfn2,Float64,false,false))
#println(@benchmark Wfn(wfn3,Float64,false,true))
println(@benchmark do_rmp2(JuWfn2))
#println(@benchmark do_direct_rmp2(JuWfn2))
println(@benchmark do_df_rmp2(JuWfn2))
println(do_rmp2(JuWfn2))
println(do_df_rmp2(JuWfn2))
#println(@benchmark do_ump2(JuWfn3))
#println(@benchmark do_ump2(JuWfn4))
