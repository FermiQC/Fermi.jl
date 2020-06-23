using Test
using PyCall
using BenchmarkTools
using Fermi.Wavefunction
using Fermi.Direct
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
mol = psi4.geometry("""
     H
     H 1 1.0
     symmetry c1
     """)
psi4.set_options(Dict("scf_type" => "pk", "d_convergence" => 14))
e, wfn = psi4.energy("hf/sto-3g", mol = mol, return_wfn = true)
JuWfn = DirectWfn(wfn)
ao_to_mo_shell(1, 1, 1, 1, JuWfn.Ca, JuWfn)
