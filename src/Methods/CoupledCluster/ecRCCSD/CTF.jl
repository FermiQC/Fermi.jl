"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function ecRCCSD{T}(Alg::CTF) where T <: AbstractFloat
    println("Generating Molecule...")
    molecule = Fermi.Geometry.Molecule()
    println("Computing Integrals...")
    aoint = Fermi.Integrals.ConventionalAOIntegrals()
    println("Calling Hartree-Fock module...")
    refwfn = Fermi.HartreeFock.RHF(molecule, aoint)

    #drop_occ = Fermi.CurrentOptions["drop_occ"]
    #drop_vir = Fermi.CurrentOptions["drop_vir"]

    #println("Transforming Integrals...")
    #moint = Fermi.Integrals.PhysRestrictedMOIntegrals{T}(refwfn.ndocc, refwfn.nvir, drop_occ, drop_vir, refwfn.C, aoint)

    #ecRCCSD{T}(refwfn, moint, Alg) 
end

"""
"""
function ecRCCSD{T}(refwfn::RHF, moint::PhysRestrictedMOIntegrals, newT1::Array{T, 2}, newT2::Array{T,4}, ecT1::Array{T,2}, ecT2::Array{T,4}, d::Array{T,2}, D::Array{T,4}, Alg::CTF) where T <: AbstractFloat

    # Print intro
    Fermi.CoupledCluster.print_header()
    @output "\n    • Computing Externally Corrected CCSD with the ecRCCSD module.\n\n"

    # Process Fock matrix, important for non HF cases
    foo = similar(moint.oo)
    foo .= moint.oo - Diagonal(moint.oo)
    fvv = similar(moint.vv)
    fvv .= moint.vv - Diagonal(moint.vv)
    fov = moint.ov
    
    # Compute Guess Energy
    Ecc = update_energy(newT1, newT2, fov, moint.oovv)
    
    @output "Initial Amplitudes Guess: MP2\n"
    @output "MP2 Energy:   {:15.10f}\n\n" Ecc+refwfn.energy

    @output "Energy from the CAS Vector:   {:15.10f}\n\n" Ecc+wfn.energy+vnuc

    # Start CC iterations
    
    cc_max_iter = Fermi.CurrentOptions["cc_max_iter"]
    cc_e_conv = Fermi.CurrentOptions["cc_e_conv"]
    cc_max_rms = Fermi.CurrentOptions["cc_max_rms"]

    @output "    Starting CC Iterations\n\n"
    @output "Iteration Options:\n"
    @output "   cc_max_iter →  {:3.0d}\n" Int(cc_max_iter)
    @output "   cc_e_conv   →  {:2.0e}\n" cc_e_conv
    @output "   cc_max_rms  →  {:2.0e}\n\n" cc_max_rms
    @output "{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}\n" "Iteration" "CC Energy" "ΔE" "Max RMS" "Time (s)"

    r1 = 1
    r2 = 1
    dE = 1
    rms = 1
    ite = 1
    T1 = similar(newT1)
    T2 = similar(newT2)

    while abs(dE) > cc_e_conv || rms > cc_max_rms
        if ite > cc_max_iter
            @output "\n⚠️  CC Equations did not converge in {:1.0d} iterations.\n" cc_max_iter
            break
        end
        t = @elapsed begin

            T1 .= newT1
            T2 .= newT2
            update_amp(T1, T2, newT1, newT2, f, V)

            # Apply external correction
            newT1 += ecT1
            newT2 += ecT2

            # Apply resolvent
            newT1 ./= d
            newT2 ./= D

            # Compute residues 
            r1 = sqrt(sum((newT1 - T1).^2))/length(T1)
            r2 = sqrt(sum((newT2 - T2).^2))/length(T2)
        end
        rms = max(r1,r2)
        oldE = Ecc
        Ecc = update_energy(newT1, newT2, f[2], V[3])
        dE = Ecc - oldE
        @output "    {:<5.0d}    {:<15.10f}    {:<12.10f}    {:<12.10f}    {:<10.5f}\n" ite Ecc dE rms t
        ite += 1
    end

    # Converged?
    if abs(dE) < cc_e_conv && rms < cc_max_rms
        @output "\n 🍾 Equations Converged!\n"
    end
    @output "\n⇒ Final ecCCSD Energy:     {:15.10f}\n" Ecc+wfn.energy+vnuc

end

