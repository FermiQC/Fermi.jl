using Fermi

# Basic options trying to be uniform

Fermi.set_num_threads()

@set {
    basis cc-pvdz
    scf_e_conv 1e-7
    scf_max_rms 1e-7
    printstyle file
    output output.dat
}

for i in 1:22

    # Read in xyz file
    molstring = read("../../xyz/S22-$(i)-dimer.xyz", String)

    # Remove first two lines from xyz
    molstring = String(split(molstring, "\n",limit=3)[end])
    Fermi.Options.set("molstring", molstring)

    # Get inthelp
    Iu = Fermi.Integrals.IntegralHelper()

    # Cold run
    wfn = @energy Iu => rhf

    # Get ingredients
    H = Iu["T"] + Iu["V"]
    F = similar(H)
    Co = wfn.orbitals.C[:, 1:wfn.ndocc]
    D = Co * Co'
    
    # Run build fock 20 times
    ttot = 0.0
    for j = 1:20
        t = @elapsed begin
            Fermi.HartreeFock.build_fock!(F, H, D, Iu)
        end
        output("@@ $j Fock Time: {}", t)
        ttot += t
    end

    # Print out time
    output("\n@@@ Average Fock Time: {}", ttot/20)

    GC.gc()
end



    

    
