"""
    Fermi.Coupled Cluster.update_energy(T1::Array{T, 2}, T2::Array{T, 4}, f::Array{T,2}, Voovv::Array{T, 4}) where T <: AbstractFloat

Compute CC energy from amplitudes and integrals.
"""
function update_energy(T1::Array{T, 2}, T2::Array{T, 4}, f::Array{T,2}, Voovv::Array{T, 4}) where { T <: AbstractFloat }

    @tensoropt (k=>x, l=>x, c=>100x, d=>100x)  begin
        CC_energy = 2.0*f[k,c]*T1[k,c]
        B[l,c,k,d] := -1.0*T1[l,c]*T1[k,d]
        B[l,c,k,d] += -1.0*T2[l,k,c,d]
        B[l,c,k,d] += 2.0*T2[k,l,c,d]
        CC_energy += B[l,c,k,d]*Voovv[k,l,c,d]
        CC_energy += 2.0*T1[l,c]*T1[k,d]*Voovv[l,k,c,d]
    end
    
    return CC_energy
end

"""
    Fermi.CoupledCluster.RCCSD.update_amp(T1::Array{T, 2}, T2::Array{T, 4}, newT1::Array{T,2}, newT2::Array{T,4}, foo::Array{T,2}, fov::Array{T,2}, fvv::Array{T,2}, moint::PhysRestrictedMOIntegrals) where T <: AbstractFloat

Update amplitudes (T1, T2) to newT1 and newT2 using CTF CCSD equations.
"""
function update_amp(T1::Array{T, 2}, T2::Array{T, 4}, newT1::Array{T,2}, newT2::Array{T,4}, foo::Array{T,2}, fov::Array{T,2}, fvv::Array{T,2}, ints::IntegralHelper, alg::A) where { T <: AbstractFloat,
                                                                                                                                                                                     A <: CCAlgorithm }

    fill!(newT1, 0.0)
    fill!(newT2, 0.0)

    # Get new amplitudes
    update_T1(T1,T2,newT1,foo,fov,fvv,ints,alg)
    update_T2(T1,T2,newT2,foo,fov,fvv,ints,alg)
end

function RCCSD{Ta}(refwfn::RHF, ints::IntegralHelper, newT1::Array{Tb, 2}, newT2::Array{Tc,4}, alg::A; ecT1 = Array{Float64}(undef,0,0), ecT2 = Array{Float64}(undef,0,0,0,0)) where { Ta <: AbstractFloat,
                                                                                                                                                                                      Tb <: AbstractFloat,
                                                                                                                                                                                      Tc <: AbstractFloat,
                                       A <: CCAlgorithm }
 
    newT1 = convert(Array{Ta},newT1)
    newT2 = convert(Array{Ta},newT2)
    # Print intro
    Fermi.CoupledCluster.print_header()
    Fermi.CoupledCluster.print_alg(alg)

    @output repeat("-",80)*"\n"
    @output "Computing and Transforming Integrals..."
    @output "\nBasis: {}\n" ints.bname["primary"]
    tint = @elapsed Fermi.CoupledCluster.compute_integrals(ints,alg)
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
    oovv = convert(Array{Ta},Fermi.CoupledCluster.compute_oovv(ints,alg))
    Ecc = update_energy(newT1, newT2, fov, oovv)
    Eguess = Ecc+refwfn.energy
    
    #@output "Initial Amplitudes Guess: MP2\n"
    @output "Guess Energy:   {:15.10f}\n\n" Ecc
    @output "Guess Total Energy:   {:15.10f}\n\n" Ecc+refwfn.energy
    
    # Start CC iterations
    
    cc_max_iter = Fermi.CurrentOptions["cc_max_iter"]
    cc_e_conv = Fermi.CurrentOptions["cc_e_conv"]
    cc_max_rms = Fermi.CurrentOptions["cc_max_rms"]
    preconv_T1 = Fermi.CurrentOptions["preconv_T1"]
    dp = Fermi.CurrentOptions["cc_damp_ratio"]
    do_diis = Fermi.CurrentOptions["diis"]
    ndiis = Fermi.CurrentOptions["ndiis"]
    do_diis ? DM_T1 = Fermi.DIIS.DIISManager{Ta,Ta}(size=6) : nothing
    do_diis ? DM_T2 = Fermi.DIIS.DIISManager{Ta,Ta}(size=6) : nothing

    @output " Dropped Occupied Orbitals →  {:3.0f}\n" Int(Fermi.CurrentOptions["drop_occ"])
    @output " Dropped Virtual Orbitals  →  {:3.0f}\n\n" Int(Fermi.CurrentOptions["drop_vir"])

    @output "Iteration Options:\n"
    @output "   cc_max_iter →  {:3.0d}\n" Int(cc_max_iter)
    @output "   cc_e_conv   →  {:2.0e}\n" cc_e_conv
    @output "   cc_max_rms  →  {:2.0e}\n\n" cc_max_rms
    @output repeat("-",80)*"\n"

    r1 = 1
    r2 = 1
    dE = 1
    rms = 1
    ite = 1
    T1 = deepcopy(newT1)
    T2 = deepcopy(newT2)

    #dfsz,nocc,nvir = size(Bov)


    @output "    Starting CC Iterations\n\n"
    preconv_T1 ? T1_time = 0 : nothing
    if preconv_T1
        @output "Preconverging T1 amplitudes\n"
        @output "Taking one T2 step\n"
        @output "{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}\n" "Iteration" "CC Energy" "ΔE" "Max RMS (T1)" "Time (s)"
        t = @elapsed begin 
            update_amp(T1, T2, newT1, newT2, foo, fov, fvv, ints, alg)

            #apply external correction
            apply_ec(newT1,zeros(Ta,size(newT2)),ecT1,ecT2)

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
        Ecc = update_energy(newT1, newT2, fov, oovv)
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
                apply_ec(newT1,zeros(Ta,size(newT2)),ecT1,ecT2)
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
                Ecc = update_energy(newT1, newT2, fov, oovv)
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

    do_diis ? DM_T1 = Fermi.DIIS.DIISManager{Ta,Ta}(size=ndiis) : nothing
    do_diis ? DM_T2 = Fermi.DIIS.DIISManager{Ta,Ta}(size=ndiis) : nothing
    if preconv_T1
        @output "Including T2 update\n"
    end

    main_time = 0
    @output "{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}\n" "Iteration" "CC Energy" "ΔE" "Max RMS" "Time (s)"

    while (abs(dE) > cc_e_conv || rms > cc_max_rms) 
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

            # Compute residues 
            r1 = sqrt(sum((newT1 .- T1).^2))/length(T1)
            r2 = sqrt(sum((newT2 .- T2).^2))/length(T2)

            if do_diis 
                e1 = (newT1 - T1)
                e2 = (newT2 - T2)
                push!(DM_T1,newT1,e1) 
                push!(DM_T2,newT2,e2) 
                if length(DM_T2) > DM_T2.max_vec
                    newT2 = Fermi.DIIS.extrapolate(DM_T2)
                    newT1 = Fermi.DIIS.extrapolate(DM_T1)
                end
            end

            newT1 .= (1-dp)*newT1 .+ dp*T1
            newT2 .= (1-dp)*newT2 .+ dp*T2
        end
        rms = max(r1,r2)
        oldE = Ecc
        Ecc = update_energy(newT1, newT2, fov, oovv)
        dE = Ecc - oldE
        main_time += t
        dE > 0 ? sign = "+" : sign = "-"
        @output "    {:<5.0d}    {:<15.10f}    {}{:>12.10f}    {:<12.10f}    {:<10.5f}\n" ite Ecc sign abs(dE) rms t
        ite += 1
    end
    @output "\nMain CCSD iterations done in {}s\n" main_time

    # Converged?
    if abs(dE) < cc_e_conv && rms < cc_max_rms 
        @output "\n 🍾 Equations Converged!\n"
    end
    @output "\n⇒ Final CCSD Energy:     {:15.10f}\n" Ecc+refwfn.energy
    @output repeat("-",80)*"\n"

    return RCCSD{Ta}(Eguess, Ecc+refwfn.energy, Fermi.MemTensor(newT1), Fermi.MemTensor(newT2))
end

function apply_ec(T1,T2,ecT1,ecT2)
    s1,s2,s3,s4 = size(ecT2)
    r1 = 1:s1
    r2 = 1:s2
    r3 = 1:s3
    r4 = 1:s4
    T1[r1,r3] .+= ecT1 
    T2[r1,r2,r3,r4] .+= ecT2
end

