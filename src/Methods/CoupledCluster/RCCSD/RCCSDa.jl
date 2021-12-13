using TensorOperations

include("RCCSDHelper.jl")

function RCCSD(moints::IntegralHelper{T,E,O}, newT1::AbstractArray{T,2}, newT2::AbstractArray{T,4}, 
                    alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
 
    # Print intro
    cc_header()
    Eref = moints.orbitals.sd_energy
    ndocc = moints.molecule.Nα
    nvir = size(moints.orbitals.C,1) - ndocc
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")

    if core ≥ ndocc
        throw(FermiException("invalid number of frozen orbitals ($core) for $ndocc doubly occupied orbitals"))
    end
    if inac ≥ nvir
        throw(FermiException("invalid number of inactive virtual orbitals ($inac) for $nvir total virtual orbitals"))
    end
    
    # Get relevant Options
    cc_max_iter = Options.get("cc_max_iter")
    cc_e_conv = Options.get("cc_e_conv")
    cc_max_rms = Options.get("cc_max_rms")
    precision_override = Options.get("precision_override")

    dp = Options.get("cc_damp_ratio")
    do_diis = Options.get("cc_diis")
    ndiis = Options.get("cc_ndiis")
    diis_start = Options.get("diis_start")
    diis_relax = Options.get("cc_diis_relax")

    diis_prec = Options.get("diis_prec")
    diis_prec = diis_prec == "double" ? Float64 : diis_prec == "single" ? Float32 : Float16
    do_diis ? DM_T1 = Fermi.DIIS.DIISManager{T, diis_prec}(size=ndiis) : nothing
    do_diis ? DM_T2 = Fermi.DIIS.DIISManager{T, diis_prec}(size=ndiis) : nothing

    output("\tDropped Occupied Orbitals →  {:3.0f}", core)
    output("\tDropped Virtual Orbitals  →  {:3.0f}", inac)
    output("\n"*repeat("-",80))
    output("\tOptions:")
    output("\tPrecision           {}", T)
    output("\tDIIS                {}", do_diis)
    output("\tDIIS Vectors        {}", (do_diis ? ndiis : 0))
    output("\tDIIS Precision      {}", diis_prec)
    output("\tDamping percentage  {}", dp)
    output("\tTBLIS               {}", Options.get("tblis"))
    output("\t\tcc_max_iter →  {:3.0d}", cc_max_iter)
    output("\t\tcc_e_conv   →  {:2.0e}", cc_e_conv)
    output("\t\tcc_max_rms  →  {:2.0e}", cc_max_rms)
    if (cc_e_conv < eps(T) || cc_max_rms < eps(T)) && !(precision_override)
        output("\t⚠ WARNING ⚠  cc_e_conv set to less than ϵ ({}) for precision {}", eps(T), T)
        output("\t             CCSD is unlikely to converge to your standards.")
        output("\t             OVERRIDING set convergence criterion.")
        output("\t             Use @set `precision_override` true to perform the computation as entered.")
        cc_e_conv = eps(T)
        cc_max_rms = eps(T)
    end

    # Initialize Loop parameters
    r1 = 1
    r2 = 1
    dE = 1
    rms = 1
    ite = 1
    T1 = deepcopy(newT1)
    T2 = deepcopy(newT2)

    output(repeat("-", 80))

    # Compute Guess Energy
    Ecc = cc_update_energy(newT1, newT2, moints, alg)
    Eguess = Ecc + Eref
    
    output("\tGuess Correlation Energy:   {:15.10f}", Ecc)
    output("\tGuess Total Energy:         {:15.10f}\n", Eguess)
    output("    Starting CC Iterations\n")

    dE = 1
    rms = 1
    ite = 1

    if do_diis
        DM_T1 = Fermi.DIIS.DIISManager{T, diis_prec}(size=ndiis)
        DM_T2 = Fermi.DIIS.DIISManager{T, diis_prec}(size=ndiis)
    end

    main_time = 0
    output("{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}    {:s5}","Iteration","CC Energy","ΔE","Max RMS","time (s)", "diis")
    while (abs(dE) > cc_e_conv || rms > cc_max_rms)
        if ite > cc_max_iter
            output("\n⚠️  CC Equations did not converge in {:1.0d} iterations.", cc_max_iter)
            break
        end
        t = @elapsed begin

            T1 .= newT1
            T2 .= newT2
            oldE = Ecc

            #debug
            update_amp!(newT1, newT2, T1, T2, moints, alg)

            # Compute residues 
            r1 = sqrt(sum((newT1 .- T1).^2)/length(T1))
            r2 = sqrt(sum((newT2 .- T2).^2)/length(T2))

            diis = false
            if do_diis
                # Update DIIS vector
                e1 = (newT1 - T1)
                e2 = (newT2 - T2)
                push!(DM_T1,newT1,e1) 
                push!(DM_T2,newT2,e2) 
                # Using DIIS at every iteration is BAD! So do it every three steps
                if ite ≥ diis_start && ite % diis_relax == 0
                    # Get new extrapolated amplitudes from DIIS
                    newT2 = Fermi.DIIS.extrapolate(DM_T2;add_res=true)
                    newT1 = Fermi.DIIS.extrapolate(DM_T1;add_res=true)
                    diis = true
                end
            end

            # Apply damping
            if dp > 0
                newT1 .= (1-dp)*newT1 .+ dp*T1
                newT2 .= (1-dp)*newT2 .+ dp*T2
            end

            Ecc = cc_update_energy(newT1, newT2, moints, alg)
        end
        rms = max(r1,r2)
        dE = Ecc - oldE
        main_time += t
        output("    {:<5.0d}    {:<15.10f}    {:>+12.10f}    {:<12.10f}    {:<10.5f}    {:s5}", ite, Ecc, dE, rms, t, diis)
        ite += 1
    end
    output("\nMain CCSD iterations done in {:5.5f} s", main_time)
    output("Average time per iteration {:5.5f}", main_time/(ite-1))

    # Converged?
    conv = false
    if abs(dE) < cc_e_conv && rms < cc_max_rms 
        output("\n 🍾 Equations Converged!")
        conv = true
    end

    output(" @Final CCSD Correlation Energy:     {:15.10f}", Ecc)
    output(" @Final CCSD Energy:                 {:15.10f}", Ecc+Eref)
    output(repeat("-",80))

    return RCCSD(Eguess, Ecc+Eref, Ecc, newT1, newT2, dE, rms)
end