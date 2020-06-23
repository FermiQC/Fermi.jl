using Test
using PyCall
using Fermi.Wavefunction
using Fermi.CISingles
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
tol = 1E-14
mol = psi4.geometry("""
                    O
                    H 1 1.1
                    H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("scf_type" => "pk", "d_convergence" => 14))
e, wfn = psi4.energy("hf/sto-3g", mol = mol, return_wfn = true)
JuWfn = Wfn(wfn)
@testset "CISingles" begin
    @testset "Smoke" begin
        out = do_RCIS(JuWfn, 5, "diag")
        @test abs(do_RCIS(JuWfn, 5)[1] - 0.320236855893771) < tol
    end
end
