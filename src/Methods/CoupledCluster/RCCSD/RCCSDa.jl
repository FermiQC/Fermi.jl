function RCCSD{T}(refwfn::RHF, ints::IntegralHelper, newT1::Array{T,2}, newT2::Array{T,4}, alg::RCCSDa) where T <: AbstractFloat
 
    # Print intro
    cc_header()

    nvir = refwfn.nvir
    ndocc = refwfn.ndocc
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")

    if core ≥ ndocc
        throw(InvalidFermiOption("invalid number of frozen orbitals ($core) for $ndocc doubly occupied orbitals"))
    end
    if inac ≥ nvir
        throw(InvalidFermiOption("invalid number of inactive virtual orbitals ($inac) for $nvir total virtual orbitals"))
    end

    # Build range for MP2 amplitudes
    # Occupied range
    o = (1+core):ndocc
    # Virtual range
    v = (ndocc+1):(ndocc + nvir - inac)

    # Collect RHF information
    orbs = refwfn.orbitals
    C = T.(orbs.C)
    ϵo = T.(orbs.eps[o])
    ϵv = T.(orbs.eps[v])
    o_size = length(ϵo)
    v_size = length(ϵv)

    if Options.get("df")
        nothing
    else
        t = @elapsed begin
        AOERI = ints["ERI"]
        Co = C[:,o]
        Cv = C[:,v]
        @tensoropt (μ=>100, ν=>100, ρ=>100, σ=>100, i=>10, j=>10, k=>10, l=>10, a=>80, b=>80, c=>80, d=>80) begin 
            OOOO[i,j,k,l] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Co[σ, l]
            OOOV[i,j,k,a] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Cv[σ, a]
            OVOV[i,a,j,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Co[ρ, j]*Cv[σ, b]
            OOVV[i,j,a,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Cv[ρ, a]*Cv[σ, b]
            OVVV[i,a,b,c] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Cv[ρ, b]*Cv[σ, c]
            OVVV[a,b,c,d] :=  AOERI[μ, ν, ρ, σ]*Cv[μ, a]*Cv[ν, b]*Cv[ρ, c]*Cv[σ, d]
        end
        AOERI = nothing
        delete!(ints, "ERI")



    tint = @elapsed Fermi.CoupledCluster.compute_integrals(ints,alg)
    delete!(ints.cache,"B") #clear out AO basis integrals for scratch space
    @output " done in {} s\n" tint
    for key in keys(ints.cache)
        ints[key] = convert(Array{Ta},ints[key])
    end

    # Process Fock matrix, important for non HF cases
    foo = convert(Array{Ta},ints["FOO"] - Diagonal(ints["FOO"]))
    fvv = convert(Array{Ta},ints["FVV"] - Diagonal(ints["FVV"]))
    fov = convert(Array{Ta},ints["FOV"])

    d = Ta[i - a for i = diag(ints["FOO"]), a = diag(ints["FVV"])]
    D = Ta[i + j - a - b for i = diag(ints["FOO"]), j = diag(ints["FOO"]), a = diag(ints["FVV"]), b = diag(ints["FVV"])]

    # Compute Guess Energy
    #oovv = convert(Array{Ta},Fermi.CoupledCluster.compute_oovv(ints,alg))
    Ecc = update_energy(newT1, newT2, fov, ints, alg)
    Eguess = Ecc+refwfn.energy
    
    #@output "Initial Amplitudes Guess: MP2\n"
    @output "\tGuess Energy:   {:15.10f}\n\n" Ecc
    @output "\tGuess Total Energy:   {:15.10f}\n\n" Ecc+refwfn.energy
    
    # Start CC iterations
    
    cc_max_iter = Fermi.CurrentOptions["cc_max_iter"]
    cc_e_conv = Fermi.CurrentOptions["cc_e_conv"]
    cc_max_rms = Fermi.CurrentOptions["cc_max_rms"]
    precision_override = Fermi.CurrentOptions["precision_override"]

    preconv_T1 = Fermi.CurrentOptions["preconv_T1"]
    dp = Fermi.CurrentOptions["cc_damp_ratio"]
    do_diis = Fermi.CurrentOptions["cc_diis"]
    ndiis = Fermi.CurrentOptions["cc_ndiis"]
    diis_start = Fermi.CurrentOptions["diis_start"]
    diis_relax = Fermi.CurrentOptions["cc_diis_relax"]
    prec_selector = Dict{String,Any}(
                                     "half" => Float16,
                                     "single" => Float32,
                                     "double" => Float64
                                    )
    diis_prec = prec_selector[Fermi.CurrentOptions["diis_prec"]]
    do_diis ? DM_T1 = Fermi.DIIS.DIISManager{Ta,diis_prec}(size=ndiis) : nothing
    do_diis ? DM_T2 = Fermi.DIIS.DIISManager{Ta,diis_prec}(size=ndiis) : nothing

    @output "\tDropped Occupied Orbitals →  {:3.0f}\n" Int(Fermi.CurrentOptions["drop_occ"])
    @output "\tDropped Virtual Orbitals  →  {:3.0f}\n\n" Int(Fermi.CurrentOptions["drop_vir"])
    @output "\n"*repeat("-",80)*"\n"
    @output "\tOptions:\n"
    @output "\tPrecision? {}\n" Ta
    @output "\tDIIS? {}\n" do_diis
    @output "\tDIIS Vectors? {}\n" (do_diis ? ndiis : 0)
    @output "\tDIIS Precision? {}\n" diis_prec
    @output "\tPreconverge T1 amplitudes? {}\n" preconv_T1
    @output "\tDamping percentage ? {}\n" dp
    @output "\t\tcc_max_iter →  {:3.0d}\n" Int(cc_max_iter)
    @output "\t\tcc_e_conv   →  {:2.0e}\n" cc_e_conv
    @output "\t\tcc_max_rms  →  {:2.0e}\n" cc_max_rms
    if (cc_e_conv < eps(Ta)) && !(precision_override)
        @output "\t⚠ WARNING⚠   cc_e_conv set to less than ϵ ({}) for precision {}\n" eps(Ta) Ta
        @output "\t             CCSD is unlikely to converge to your standards.\n"
        @output "\t             OVERRIDING set convergence criterion.\n"
        @output "\t             Use @set precision_override true to perform the computation as entered.\n"
        cc_e_conv = eps(Ta)
    end
    if (cc_max_rms < eps(Ta)) && !(precision_override)
        @output "\t⚠ WARNING⚠   cc_max_rms set to less than ϵ ({}) for precision {}\n" eps(Ta) Ta
        @output "\t             CCSD is unlikely to converge to your standards.\n"
        @output "\t             OVERRIDING set convergence criterion.\n"
        @output "\t             Use @set precision_override true to perform the computation as entered.\n"
        cc_max_rms = eps(Ta)
    end

    r1 = 1
    r2 = 1
    dE = 1
    rms = 1
    ite = 1
    T1 = deepcopy(newT1)
    T2 = deepcopy(newT2)

    #dfsz,nocc,nvir = size(Bov)


    @output repeat("-",80)*"\n"
    @output "    Starting CC Iterations\n\n"
    preconv_T1 ? T1_time = 0 : nothing
    dummyT2 = zeros(Ta,size(newT2))
    if preconv_T1
        @output "Preconverging T1 amplitudes\n"
        @output "Taking one T2 step\n"
        @output "{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}\n" "Iteration" "CC Energy" "ΔE" "Max RMS (T1)" "Time (s)"
        t = @elapsed begin 
            update_amp(T1, T2, newT1, newT2, foo, fov, fvv, ints, alg)

            #apply external correction
            apply_ec(newT1,dummyT2,ecT1,ecT2)

            # Apply resolvent
            newT1 ./= d
            newT2 ./= D

            # Compute residues 
            r1 = sqrt(sum((newT1 .- T1).^2))/length(T1)
            r2 = sqrt(sum((newT2 .- T2).^2))/length(T2)

            if do_diis 
                e1 = (newT1 - T1)
                e2 = (newT2 - T2)
                push!(DM_T1,newT1,e1) 
                push!(DM_T2,newT2,e2) 
            end

            newT1 .= (1-dp)*newT1 .+ dp*T1
            newT2 .= (1-dp)*newT2 .+ dp*T2
        end
        T1_time += t

        rms = max(r1,r2)
        oldE = Ecc
        Ecc = update_energy(newT1, newT2, fov, ints, alg)
        dE = Ecc - oldE
        @output "    {:<5}    {:<15.10f}    {:<12.10f}    {:<12.10f}    {:<10.5f}\n" "pre" Ecc dE rms t

        while abs(dE) > cc_e_conv || rms > cc_max_rms
            if ite > cc_max_iter
                @output "\n⚠️  CC Equations did not converge in {:1.0d} iterations.\n" cc_max_iter
                break
            end
            t = @elapsed begin
                T1 .= newT1
                T2 .= newT2
                update_T1(T1,T2,newT1,foo,fov,fvv,ints,alg)
                apply_ec(newT1,dummyT2,ecT1,ecT2)
                newT1 ./= d
                if do_diis 
                    e1 = newT1 - T1
                    push!(DM_T1,newT1,e1) 
                    newT1 = Fermi.DIIS.extrapolate(DM_T1)
                end

                # Compute residues 
                r1 = sqrt(sum((newT1 .- T1).^2))/length(T1)

                newT1 .= (1-dp)*newT1 .+ dp*T1

                rms = r1
                oldE = Ecc
                Ecc = update_energy(newT1, newT2, fov, ints, alg)
                dE = Ecc - oldE
            end
            T1_time += t
            dE > 0 ? sign = "+" : sign = "-"
            @output "    {:<5.0d}    {:<15.10f}    {}{:>12.10f}    {:<12.10f}    {:<10.5f}\n" ite Ecc sign abs(dE) rms t
            ite += 1
        end
        @output "\nT1 pre-convergence took {}s\n" T1_time
    end

    dE = 1
    rms = 1
    ite = 1

    do_diis ? DM_T1 = Fermi.DIIS.DIISManager{Ta,diis_prec}(size=ndiis) : nothing
    do_diis ? DM_T2 = Fermi.DIIS.DIISManager{Ta,diis_prec}(size=ndiis) : nothing
    if preconv_T1
        @output "Including T2 update\n"
    end

    main_time = 0
    @output "{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}\n" "Iteration" "CC Energy" "ΔE" "Max RMS" "Time (s)"
    relax = 1
            Ecc = update_energy(newT1, newT2, fov, ints, alg)
            @output "{}\n" Ecc

    while (abs(dE) > cc_e_conv || rms > cc_max_rms) || ite < 10
        if ite > cc_max_iter
            @output "\n⚠️  CC Equations did not converge in {:1.0d} iterations.\n" cc_max_iter
            break
        end
        t = @elapsed begin

            T1 .= newT1
            T2 .= newT2
            update_amp(T1, T2, newT1, newT2, foo, fov, fvv, ints, alg)

            apply_ec(newT1,newT2,ecT1,ecT2)

            # Apply resolvent
            newT1 ./= d
            newT2 ./= D
            oldE = Ecc
            Ecc = update_energy(newT1, newT2, fov, ints, alg)

            # Compute residues 
            r1 = sqrt(sum((newT1 .- T1).^2))/length(T1)
            r2 = sqrt(sum((newT2 .- T2).^2))/length(T2)

            if do_diis && ite > 2
                relax -= 1
                e1 = (newT1 - T1)
                e2 = (newT2 - T2)
                push!(DM_T1,newT1,e1) 
                push!(DM_T2,newT2,e2) 
                if relax == 0 
                    newT2 = Fermi.DIIS.extrapolate(DM_T2;add_res=true)
                    newT1 = Fermi.DIIS.extrapolate(DM_T1;add_res=true)
                    relax = diis_relax
                end
            end

            newT1 .= (1-dp)*newT1 .+ dp*T1
            newT2 .= (1-dp)*newT2 .+ dp*T2
        end
        rms = max(r1,r2)
        Ecc = update_energy(newT1, newT2, fov, ints, alg)
        dE = Ecc - oldE
        main_time += t
        dE > 0 ? sign = "+" : sign = "-"
        @output "    {:<5.0d}    {:<15.10f}    {}{:>12.10f}    {:<12.10f}    {:<10.5f}\n" ite Ecc sign abs(dE) rms t
        ite += 1
    end
    @output "\nMain CCSD iterations done in {}s\n" main_time

    # Converged?
    conv = false
    if abs(dE) < cc_e_conv && rms < cc_max_rms 
        @output "\n 🍾 Equations Converged!\n"
        conv = true
    end
    @output "\n⇒ Final CCSD Energy:     {:15.10f}\n" Ecc+refwfn.energy
    @output repeat("-",80)*"\n"

    return RCCSD{Ta}(Eguess, Ecc+refwfn.energy, Fermi.MemTensor(newT1), Fermi.MemTensor(newT2), conv)
end