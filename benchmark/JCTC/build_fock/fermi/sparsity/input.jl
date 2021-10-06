using Fermi
using PrettyTables

# Basic options trying to be uniform

Fermi.set_num_threads()

@set {
    basis cc-pvdz
    scf_e_conv 1e-7
    scf_max_rms 1e-7
    printstyle file
    output output.dat
}

S = zeros(22)

for i in 1:22

    # Read in xyz file
    molstring = read("../../../xyz/S22-$(i)-dimer.xyz", String)

    # Remove first two lines from xyz
    molstring = String(split(molstring, "\n",limit=3)[end])
    Fermi.Options.set("molstring", molstring)
    output("Molecule: S$i")

    # Get inthelp
    Iu = Fermi.Integrals.IntegralHelper()

    eri = Iu["ERI"]

    N::Int128 = Iu.orbitals.basisset.nbas
    output("Number of basis: {}", N)
    output("Unique ERI elements: {}", length(eri.data))


    S[i] = 100* (length(eri.data) / N^4)
    output("Sparsity: {}", S[i])

    GC.gc()
end

@pt :header = ["Molecule", "Sparsity"] [collect(1:22) S]
