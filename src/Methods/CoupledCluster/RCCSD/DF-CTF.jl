using LinearAlgebra
using Fermi.DIIS

"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function RCCSD{T}(guess::RCCSD{Tb},Alg::DFCTF) where { T <: AbstractFloat,
                                                    Tb <: AbstractFloat }
    molecule = Fermi.Geometry.Molecule()
    aoint = Fermi.Integrals.DFAOIntegrals(molecule)
    refwfn = Fermi.HartreeFock.RHF(molecule, aoint, Fermi.HartreeFock.DFRHF())

    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    @output "Transforming Integrals..."
    tint = @elapsed moint = Fermi.Integrals.DFMOIntegrals{T}(refwfn.ndocc, refwfn.nvir, drop_occ, drop_vir, refwfn.C, aoint)
    @output " done in {} s" tint
    RCCSD{T}(refwfn, guess, moint, Alg) 
end

"""
    Fermi.CoupledCluster.RCCSD{T}(refwfn::RHF, moint::PhysRestrictedMOIntegrals, Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF). Precision (T), reference wavefunction (refwfn)
and molecular orbital integrals (moint) must be passed.
"""
function RCCSD{T}(refwfn::RHF, moint::PhysRestrictedMOIntegrals, Alg::DFCTF) where T <: AbstractFloat
    d = [i - a for i = diag(moint.oo), a = diag(moint.vv)]
    D = [i + j - a - b for i = diag(moint.oo), j = diag(moint.oo), a = diag(moint.vv), b = diag(moint.vv)]
    newT1 = moint.ov./d
    Bov = moint.Bov
    @tensor oovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
    #oovv = permutedims(zeros(size(D)),(1,3,2,4))
    #Fermi.contract!(oovv,aoint.Bov,aoint.Bov,"pqrs","Qpq","Qrs")
    #oovv = permutedims(oovv,(1,3,2,4))
    newT2 = oovv./D
    RCCSD{T}(refwfn, moint, newT1, newT2, Alg)
end

function RCCSD{T}(refwfn::RHF, guess::RCCSD{Tb}, moint::Fermi.Integrals.DFMOIntegrals, Alg::DFCTF) where { T <: AbstractFloat,
                                                                                                    Tb <: AbstractFloat }
    d = [i - a for i = diag(moint.oo), a = diag(moint.vv)]
    D = [i + j - a - b for i = diag(moint.oo), j = diag(moint.oo), a = diag(moint.vv), b = diag(moint.vv)]
    newT1 = moint.ov./d
    Bov = moint.Bov
    @tensor oovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
    println(size(oovv))
    #oovv = permutedims(zeros(size(D)),(1,3,2,4))
    #Fermi.contract!(oovv,moint.Bov,moint.Bov,"pqrs","Qpq","Qrs")
    #oovv = permutedims(oovv,(1,3,2,4))
    newT2 = oovv ./ D
    @output "Here1\n"
    RCCSD{T}(refwfn, moint, newT1, newT2, Alg)
end

"""
    Fermi.CoupledCluster.RCCSD{T}(refwfn::RHF, moint::PhysRestrictedMOIntegrals, Alg::CTF)

Base function for CTF RCCSD.
"""
function RCCSD{T}(refwfn::RHF, moint::Fermi.Integrals.DFMOIntegrals, newT1::Array{T, 2}, newT2::Array{T,4}, Alg::DFCTF) where T <: AbstractFloat

    # Print intro
    Fermi.CoupledCluster.print_header()
    @output "\n    • Computing CCSD with the DF-CFT algorithm .\n\n"
    d = [i - a for i = diag(moint.oo), a = diag(moint.vv)]
    D = [i + j - a - b for i = diag(moint.oo), j = diag(moint.oo), a = diag(moint.vv), b = diag(moint.vv)]

    # Process Fock matrix, important for non HF cases
    foo = similar(moint.oo)
    foo .= moint.oo - Diagonal(moint.oo)
    fvv = similar(moint.vv)
    fvv .= moint.vv - Diagonal(moint.vv)
    fov = moint.ov

    # Compute Guess Energy
    Bov = moint.Bov
    @tensor oovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
    Ecc = update_energy(newT1, newT2, fov, oovv)
    Eguess = Ecc+refwfn.energy
    
    @output "Initial Amplitudes Guess: MP2\n"
    @output "MP2 Energy:   {:15.10f}\n\n" Ecc
    @output "MP2 Total Energy:   {:15.10f}\n\n" Ecc+refwfn.energy
    
    # Start CC iterations
    
    cc_max_iter = Fermi.CurrentOptions["cc_max_iter"]
    cc_e_conv = Fermi.CurrentOptions["cc_e_conv"]
    cc_max_rms = Fermi.CurrentOptions["cc_max_rms"]
    preconv_T1 = Fermi.CurrentOptions["preconv_T1"]
    dp = Fermi.CurrentOptions["cc_damp_ratio"]
    do_diis = Fermi.CurrentOptions["diis"]
    do_diis ? DM_T1 = Fermi.DIIS.DIISManager{Float64,Float64}(size=8) : nothing
    do_diis ? DM_T2 = Fermi.DIIS.DIISManager{Float64,Float64}(size=8) : nothing


    @output "    Starting CC Iterations\n\n"
    @output "Iteration Options:\n"
    @output "   cc_max_iter →  {:3.0d}\n" Int(cc_max_iter)
    @output "   cc_e_conv   →  {:2.0e}\n" cc_e_conv
    @output "   cc_max_rms  →  {:2.0e}\n\n" cc_max_rms

    r1 = 1
    r2 = 1
    dE = 1
    rms = 1
    ite = 1
    T1 = deepcopy(newT1)
    T2 = deepcopy(newT2)

    dfsz,nocc,nvir = size(moint.Bov)
    T2_inter = zeros(T,dfsz,nocc,nocc,nvir,nvir)


    preconv_T1 ? T1_time = 0 : nothing
    if preconv_T1
        @output "Preconverging T1 amplitudes\n"
        @output "Taking one T2 step\n"
        @output "{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}\n" "Iteration" "CC Energy" "ΔE" "Max RMS (T1)" "Time (s)"
        t = @elapsed begin 
            update_amp(T1, T2, newT1, newT2, foo, fov, fvv, moint,T2_inter)

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
                #newT1 = Fermi.DIIS.extrapolate(DM_T1)
                #newT2 = Fermi.DIIS.extrapolate(DM_T2)
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
                update_T1(T1,T2,newT1,foo,fov,fvv,moint)
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
            @output "    {:<5.0d}    {:<15.10f}    {:<12.10f}    {:<12.10f}    {:<10.5f}\n" ite Ecc dE rms t
            ite += 1
        end
        @output "\nT1 pre-convergence took {}s\n" T1_time
    end

    dE = 1
    rms = 1
    ite = 1

    do_diis ? DM_T1 = Fermi.DIIS.DIISManager{Float64,Float64}(size=6) : nothing
    do_diis ? DM_T2 = Fermi.DIIS.DIISManager{Float64,Float64}(size=6) : nothing
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
            update_amp(T1, T2, newT1, newT2, foo, fov, fvv, moint,T2_inter)

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
        @output "    {:<5.0d}    {:<15.10f}    {:<12.10f}    {:<12.10f}    {:<10.5f}\n" ite Ecc dE rms t
        ite += 1
    end
    @output "\nMain CCSD iterations done in {}s\n" main_time

    # Converged?
    if abs(dE) < cc_e_conv && rms < cc_max_rms 
        @output "\n 🍾 Equations Converged!\n"
    end
    @output "\n⇒ Final CCSD Energy:     {:15.10f}\n" Ecc+refwfn.energy

    return RCCSD{T}(Eguess, Ecc+refwfn.energy, Fermi.MemTensor(newT1), Fermi.MemTensor(newT2))
end

"""
    Fermi.Coupled Cluster.update_energy(T1::Array{T, 2}, T2::Array{T, 4}, f::Array{T,2}, Voovv::Array{T, 4}) where T <: AbstractFloat

Compute CC energy from amplitudes and integrals.
"""
function update_energy(T1::Array{T, 2}, T2::Array{T, 4}, f::Array{T,2}, Voovv::Array{T, 4}) where T <: AbstractFloat

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
function update_amp(T1::Array{T, 2}, T2::Array{T, 4}, newT1::Array{T,2}, newT2::Array{T,4}, foo::Array{T,2}, fov::Array{T,2}, fvv::Array{T,2}, moint::Fermi.Integrals.DFMOIntegrals, T2_inter) where T <: AbstractFloat

    #Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = moint.oooo, moint.ooov, moint.oovv, moint.ovov, moint.ovvv, moint.vvvv

    fill!(newT1, 0.0)
    fill!(newT2, 0.0)
    fill!(T2_inter, 0.0)

    # Get new amplitudes
    update_T1(T1,T2,newT1,foo,fov,fvv,moint)
    @output "T1 done\n"
    update_T2(T1,T2,newT2,foo,fov,fvv,moint,T2_inter)
end

function update_T1(T1::Array{T,2}, T2::Array{T,4}, newT1::Array{T,2}, foo, fov, fvv, moint::Fermi.Integrals.DFMOIntegrals) where T <: AbstractFloat
    #Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = moint.oooo, moint.ooov, moint.oovv, moint.ovov, moint.ovvv, moint.vvvv
    Bov = moint.Bov
    Bvo = moint.Bvo
    Boo = moint.Boo
    Bvv = moint.Bvv
    @tensor begin
        Voooo[i,j,k,l] := Boo[Q,i,k]*Boo[Q,j,l]
        Vooov[i,j,k,a] := Boo[Q,i,k]*Bov[Q,j,a]
        Voovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
        Vovov[i,a,j,b] := Boo[Q,i,j]*Bvv[Q,a,b]
        #Vovvv[i,a,b,c] := Bov[Q,i,b]*Bvv[Q,a,c]
        #Vvvvv[a,b,c,d] := Bvv[Q,a,c]*Bvv[Q,b,d]
    end
    #@tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
    #    newT1[i,a] += fov[i,a]
    #    newT1[i,a] -= foo[i,k]*T1[k,a]
    #    newT1[i,a] += fvv[c,a]*T1[i,c]
    #    newT1[i,a] -= fov[k,c]*T1[i,c]*T1[k,a]
    #    newT1[i,a] += 2.0*fov[k,c]*T2[i,k,a,c]
    #    newT1[i,a] -= fov[k,c]*T2[k,i,a,c]
    #    newT1[i,a] -= T1[k,c]*Vovov[i,c,k,a]
    #    newT1[i,a] += 2.0*T1[k,c]*Voovv[k,i,c,a]
    #    #newT1[i,a] -= T2[k,i,c,d]*Vovvv[k,a,d,c]
    #    newT1[i,a] -= T2[k,i,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
    #    #newT1[i,a] += 2.0*T2[i,k,c,d]*Vovvv[k,a,d,c]
    #    newT1[i,a] += 2.0*T2[i,k,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
    #    newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,l,i,c]
    #    newT1[i,a] += T2[l,k,a,c]*Vooov[k,l,i,c]
    #    newT1[i,a] += -2.0*T1[k,c]*T1[l,a]*Vooov[l,k,i,c]
    #    #newT1[i,a] -= T1[k,c]*T1[i,d]*Vovvv[k,a,d,c]
    #    newT1[i,a] -= T1[k,c]*T1[i,d]*Bov[Q,k,d]*Bvv[Q,a,c]
    #    #newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Vovvv[k,a,c,d]
    #    newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Bov[Q,k,c]*Bvv[Q,a,d]
    #    newT1[i,a] += T1[k,c]*T1[l,a]*Vooov[k,l,i,c]
    #    newT1[i,a] += -2.0*T1[k,c]*T2[i,l,a,d]*Voovv[l,k,c,d]
    #    newT1[i,a] += -2.0*T1[k,c]*T2[l,i,a,d]*Voovv[k,l,c,d]
    #    newT1[i,a] += T1[k,c]*T2[l,i,a,d]*Voovv[l,k,c,d]
    #    newT1[i,a] += -2.0*T1[i,c]*T2[l,k,a,d]*Voovv[l,k,c,d]
    #    newT1[i,a] += T1[i,c]*T2[l,k,a,d]*Voovv[k,l,c,d]
    #    newT1[i,a] += -2.0*T1[l,a]*T2[i,k,d,c]*Voovv[k,l,c,d]
    #    newT1[i,a] += T1[l,a]*T2[i,k,c,d]*Voovv[k,l,c,d]
    #    newT1[i,a] += T1[k,c]*T1[i,d]*T1[l,a]*Voovv[l,k,c,d]
    #    newT1[i,a] += -2.0*T1[k,c]*T1[i,d]*T1[l,a]*Voovv[k,l,c,d]
    #    newT1[i,a] += 4.0*T1[k,c]*T2[i,l,a,d]*Voovv[k,l,c,d]
    #end
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x,Q=>50x) begin
        newT1[i,a] += fov[i,a]
        newT1[i,a] -= foo[i,k]*T1[k,a]
        newT1[i,a] += fvv[c,a]*T1[i,c]
        newT1[i,a] -= fov[k,c]*T1[i,c]*T1[k,a]
        newT1[i,a] += 2.0*fov[k,c]*T2[i,k,a,c]
        newT1[i,a] -= fov[k,c]*T2[k,i,a,c]
        newT1[i,a] -= T1[k,c]*Vovov[i,c,k,a]
        newT1[i,a] += 2.0*T1[k,c]*Voovv[k,i,c,a]
        #newT1[i,a] -= T2[k,i,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] -= T2[k,i,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        #newT1[i,a] += 2.0*T2[i,k,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += T2[l,k,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += -2.0*T1[k,c]*T1[l,a]*Vooov[l,k,i,c]
        #newT1[i,a] -= T1[k,c]*T1[i,d]*Vovvv[k,a,d,c]
        newT1[i,a] -= T1[k,c]*T1[i,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        #newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Vovvv[k,a,c,d]
        newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Bov[Q,k,c]*Bvv[Q,a,d]
        newT1[i,a] += T1[k,c]*T1[l,a]*Vooov[k,l,i,c]
        newT1[i,a] += -2.0*T1[k,c]*T2[i,l,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[k,c]*T2[l,i,a,d]*Voovv[k,l,c,d]
        newT1[i,a] += T1[k,c]*T2[l,i,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[i,c]*T2[l,k,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += T1[i,c]*T2[l,k,a,d]*Voovv[k,l,c,d]
        newT1[i,a] += -2.0*T1[l,a]*T2[i,k,d,c]*Voovv[k,l,c,d]
        newT1[i,a] += T1[l,a]*T2[i,k,c,d]*Voovv[k,l,c,d]
        newT1[i,a] += T1[k,c]*T1[i,d]*T1[l,a]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[k,c]*T1[i,d]*T1[l,a]*Voovv[k,l,c,d]
        newT1[i,a] += 4.0*T1[k,c]*T2[i,l,a,d]*Voovv[k,l,c,d]
    end
end

function update_T2(T1::Array{T,2},T2::Array{T,4},newT2::Array{T,4},foo,fov,fvv,moint::Fermi.Integrals.DFMOIntegrals, T2_inter) where T <: AbstractFloat
    #Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = moint.oooo, moint.ooov, moint.oovv, moint.ovov, moint.ovvv, moint.vvvv
    Bov = moint.Bov
    Bvo = moint.Bvo
    Boo = moint.Boo
    Bvv = moint.Bvv
    @tensor begin
        Voooo[i,j,k,l] := Boo[Q,i,k]*Boo[Q,j,l]
        Vooov[i,j,k,a] := Boo[Q,i,k]*Bov[Q,j,a]
        Voovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
        Vovov[i,a,j,b] := Boo[Q,i,j]*Bvv[Q,a,b]
        #Vovvv[i,a,b,c] := Bov[Q,i,b]*Bvv[Q,a,c]
        #Vvvvv[a,b,c,d] := Bvv[Q,a,c]*Bvv[Q,b,d]
    end
    dfsz,nocc,nvir = size(Bov)
    #T2_inter = zeros(dfsz,nocc,nocc,nvir,nvir)
    @tensor _T[i,j,c,d] := T1[i,c]*T1[j,d] + T2[i,j,c,d]
    #Fermi.contract!(T2_inter,_T,Bvv,"Qijad","ijcd","Qca")
    #Fermi.contract!(newT2,T2_inter,Bvv,"ijab","Qijad","Qdb")
    T2_inter = zeros(dfsz,nvir,nvir)
    T2_slice = zeros(nvir,nvir)
    Tslice = zeros(nvir,nvir)
    for i=1:nocc, j=1:nocc
        T2_inter .= 0
        T2_slice .= 0
        Tslice .= _T[i,j,:,:]
        Fermi.contract!(T2_inter,Tslice,Bvv,"Qad","cd","Qca")
        Fermi.contract!(T2_slice,T2_inter,Bvv,"ab","Qad","Qdb")
        newT2[i,j,:,:] += T2_slice
    end


    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>50x) begin
        newT2[i,j,a,b] += Voovv[i,j,a,b]
        #newT2[i,j,a,b] += Bvv[Q,c,a]*T1[i,c]*T1[j,d]*Bvv[Q,d,b]
        #newT2[i,j,a,b] += Bvv[Q,c,a]*T2[i,j,c,d]*Bvv[Q,d,b]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] += T2[k,l,a,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,a]*Bov[Q,k,c]*Bvv[Q,b,d]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,b]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT2[i,j,a,b] += T1[i,c]*T1[k,a]*T1[l,b]*Vooov[l,k,j,c]
        newT2[i,j,a,b] += T1[j,c]*T1[k,a]*T1[l,b]*Vooov[k,l,i,c]
        newT2[i,j,a,b] += T2[k,l,a,c]*T2[i,j,d,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[l,j,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[l,k,a,c]*T2[i,j,d,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,d,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T2[i,k,a,c]*T2[l,j,b,d]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[j,l,b,d]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[k,i,a,c]*T2[j,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[i,j,a,c]*T2[l,k,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,j,a,c]*T2[k,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[k,j,a,c]*T2[i,l,d,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += 4.0*T2[i,k,a,c]*T2[j,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[i,j,d,c]*T2[l,k,a,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T1[k,a]*T1[l,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T2[l,k,a,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*T2[i,j,d,c]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] := -1.0*foo[i,k]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += fvv[c,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[k,b]*Vooov[j,i,k,a]
        P_OoVv[i,j,a,b] += T1[j,c]*Bov[Q,i,a]*Bvv[Q,c,b]
        P_OoVv[i,j,a,b] += -1.0*fov[k,c]*T1[i,c]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += -1.0*fov[k,c]*T1[k,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,i,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,a]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,b]*Vovov[j,c,k,a]
        P_OoVv[i,j,a,b] += 2.0*T2[i,k,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[i,k,a,c]*Vovov[j,c,k,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,j,a,c]*Vovov[i,c,k,b]
        P_OoVv[i,j,a,b] += -2.0*T1[l,b]*T2[i,k,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[k,i,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,d,b]*Bov[Q,k,c]*Bvv[Q,a,d]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[k,i,a,d]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += (-1.0*T1[j,c]*T2[i,k,a,d] + 2.0*T1[k,c]*T2[i,j,a,d])*Bov[Q,k,c]*Bvv[Q,b,d]
        P_OoVv[i,j,a,b] += T1[j,c]*T2[l,k,a,b]*Vooov[l,k,i,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[i,k,a,c]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,a]*T2[i,j,d,c]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += T1[k,a]*T2[i,l,c,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += (2.0*T1[j,c]*T2[i,k,a,d] - 1.0*T1[k,c]*T2[i,j,a,d])*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += T1[k,c]*T2[i,l,a,b]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T2[i,l,a,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += T2[j,k,c,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[i,c]*T1[k,a]*T2[j,l,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[l,b]*T2[k,j,a,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T2[i,k,d,c]*T2[l,j,a,b]*Voovv[k,l,c,d]
        
        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end
end
