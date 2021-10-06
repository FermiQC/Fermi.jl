using Fermi
using Fermi.Integrals
using MKL

# Cold run
@set printstyle none
@set df true
@energy mp2

# Basic options trying to be uniform

Fermi.set_num_threads()

@set {
    basis cc-pvdz
    rifit cc-pvtz-rifit
    jkfit cc-pvtz-jkfit
    scf_e_conv 1e-8
    scf_max_rms 1e-8
    printstyle file
    output output.dat
}

for i in 1:22

    # Read in xyz file
    molstring = read("../../xyz/S22-$(i)-dimer.xyz", String)

    # Remove first two lines from xyz
    molstring = String(split(molstring, "\n",limit=3)[end])
    Fermi.Options.set("molstring", molstring)


    # Compute SCF
    wfn = @energy rhf

    # Get AO Integral helper
    aoints = IntegralHelper(eri_type=RIFIT())

    # Compute ERI
    aoints["ERI"]

    # Get MO Integral helper
    moints = IntegralHelper(orbitals=wfn.orbitals, eri_type=RIFIT())

    # Time integral transformation
    timings = zeros(10)
    output("Timing integral transformation")
    for n = 1:10
        t = @elapsed Fermi.Integrals.compute_BOV!(moints, aoints)
        timings[n] = t
        output("Integral Transf. S$i  $n run: $t")
    end
    output("@@ Average time for integral transformation: {:15.10f}", sum(timings)/10)
    
    # Time energy computation
    timings = zeros(10)
    output("Timing DF-MP2 energy computation")
    for n = 1:10
        t = @elapsed @energy moints => mp2
        timings[n] = t
        output("DF-MP2 Energy: S$i  $n run: $t")
    end
    output("@@ Average time for DF-MP2 energy: {:15.10f}", sum(timings)/10)

    GC.gc()
end



    

    
