using Fermi.Integrals

using Fermi.ConfigurationInteraction.DetOperations

function ecRCCSD(alg::ecRCCSDa)

    aoints = IntegralHelper{Float64}()
    rhf = Fermi.HartreeFock.RHF(aoints)

    ciwfn = Fermi.ConfigurationInteraction.RFCI(aoints, rhf)

    if ciwfn.dets === nothing
        throw(ArgumentError("CI Wave function did not return determinants. Please select a different algorithm."))
    end

    # Process CI wave function
    refdet, Casdata = process_cas(ciwfn) 

    # Transform Integrals
    if typeof(aoints.eri_type) === JKFIT || Options.get("precision") == "single"
        aoints = IntegralHelper(eri_type=RIFIT())
    end
    moints = IntegralHelper(orbitals=rhf.orbitals)
    Fermi.Integrals.compute!(moints, aoints, "F")
    Fermi.Integrals.compute!(moints, aoints, "OOOO")
    Fermi.Integrals.compute!(moints, aoints, "OOOV")
    Fermi.Integrals.compute!(moints, aoints, "OOVV")
    Fermi.Integrals.compute!(moints, aoints, "OVOV")
    Fermi.Integrals.compute!(moints, aoints, "OVVV")
    Fermi.Integrals.compute!(moints, aoints, "VVVV")

    drop_occ = Fermi.Options.get("drop_occ")

    # Perform cluster cluster decomposition
    output("   • CAS Decomposition started:")

    # NOTE: The cas decomposition routine was originally coded for integrals in Physicists' notation
    # Our latest Fermi version uses chemists' notation. Hence, here before calling the cas_decomposition function
    # We do a permutation to input physiscist integrals... It is not an ideal situation, but such is life.
    p(X) = permutedims(X, (1,3,2,4))
    t = @elapsed T1, T2, ecT1, ecT2 = cas_decomposition(refdet, Casdata, rhf.ndocc, drop_occ, moints["Fia"], p(moints["OVOV"]), p(moints["OVVV"]), p(moints["OOOV"]))
    println("aya")
    output("Finished in {:5.5} seconds.\n", t)
    Casdata = nothing

    ecRCCSD(moints, T1, T2, ecT1, ecT2, alg)
end

function ecRCCSD(moints::IntegralHelper{T,E,O}, T1, T2, ecT1, ecT2, alg::ecRCCSDa) where {T<:AbstractFloat,E<:AbstractERI,O<:AbstractRestrictedOrbitals}

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
    newT1 = deepcopy(T1)
    newT2 = deepcopy(T2)

    output(repeat("-", 80))

    # Compute Guess Energy
    Ecc = cc_update_energy(newT1, newT2, moints, RCCSDa())
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

            update_amp!(newT1, newT2, T1, T2, moints, ecT1, ecT2, alg)

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

            Ecc = cc_update_energy(newT1, newT2, moints, RCCSDa())
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

    output(" @Final ecCCSD Correlation Energy:     {:15.10f}", Ecc)
    output(" @Final ecCCSD Energy:                 {:15.10f}", Ecc+Eref)
    output(repeat("-",80))

    if Options.get("return_ints")
        return RCCSD(Eguess, Ecc, Ecc+Eref, newT1, newT2, dE, rms), moints
    else
        return RCCSD(Eguess, Ecc, Ecc+Eref, newT1, newT2, dE, rms)
    end
end

function update_amp!(newT1::AbstractArray{T,2}, newT2::AbstractArray{T,4}, T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,E,RHFOrbitals}, 
                    ecT1::AbstractArray{T, 2}, ecT2::AbstractArray{T, 4}, alg::ecRCCSDa) where {T<:AbstractFloat,E<:AbstractERI}

    # Clean the arrays
    fill!(newT1, 0.0)
    fill!(newT2, 0.0)

    # Get new amplitudes
    cc_update_T1!(newT1, T1, T2, moints, RCCSDa())
    cc_update_T2!(newT2, T1, T2, moints, RCCSDa())

    # Add external correction
    newT1 .+= ecT1
    newT2 .+= ecT2

    # Orbital energies line
    if haskey(moints.cache, "D1")
        d = moints["D1"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        frozen = Options.get("drop_occ")
        inac = Options.get("drop_vir")
        ϵo = Fd[(1+frozen:ndocc)]
        ϵv = Fd[(1+ndocc):end-inac]

        d = [ϵo[i]-ϵv[a] for i=eachindex(ϵo), a=eachindex(ϵv)]
        moints["D1"] = d
    end

    if haskey(moints.cache, "D2")
        D = moints["D2"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        frozen = Options.get("drop_occ")
        inac = Options.get("drop_vir")
        ϵo = Fd[(1+frozen:ndocc)]
        ϵv = Fd[(1+ndocc):end-inac]

        D = [ϵo[i]+ϵo[j]-ϵv[a]-ϵv[b] for i=eachindex(ϵo), j=eachindex(ϵo), a=eachindex(ϵv), b=eachindex(ϵv)]
        moints["D2"] = D
    end

    # Orbital energies line
    d, D = moints["D1"], moints["D2"]

    newT1 ./= d
    newT2 ./= D
end