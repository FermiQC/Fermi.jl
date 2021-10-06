temp = """
import psi4

# Basic options trying to be uniform

psi4.set_num_threads(10)
psi4.set_memory('90 GB')
psi4.set_options({"basis" : "cc-pvdz",
                  "scf_type" : "df",
                  "cc_type" : "df",
                  "df_basis_scf" : "cc-pvtz-jkfit",
                  "df_basis_cc" : "cc-pvtz-ri",
                  "puream" : True,
                  "TRIPLES_LOW_MEMORY" : False,
                  "e_convergence" : 1e-8,
                  "d_convergence" : 1e-8})

molecule {
[molecule]
symmetry c1
}

energy("ccsd(t)")
"""

for i = 1:22
    xyz = readlines("../../xyz/S22-$i-dimer.xyz")[3:end]
    xyz = join(xyz, "\n")

    mkdir("S$i")
    open("S$i/input.dat", "w") do io
        write(io, replace(temp, "[molecule]"=>xyz))
    end
end
