using LoopVectorization
using LinearAlgebra
using TensorOperations

function RMP2()
    aoints = IntegralHelper{Float64}()
    rhf = RHF(aoints)

    if typeof(aoints.eri_type) === JKFIT || Options.get("precision") == "single"
        aoints = IntegralHelper(eri_type=RIFIT())
    end
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RMP2(moints, aoints)
end

function RMP2(aoints::IntegralHelper{Float64,E,AtomicOrbitals}) where E <: AbstractERI
    rhf = RHF(aoints)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RMP2(moints, aoints)
end

function RMP2(O::AbstractRestrictedOrbitals)
    moints = IntegralHelper(orbitals=O)
    aoints = IntegralHelper(eri_type=moints.eri_type)
    RMP2(moints, aoints)
end

function RMP2(rhf::RHF)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    aoints = IntegralHelper(eri_type=moints.eri_type)
    RMP2(moints, aoints)
end

function RMP2(M::Molecule)
    aoints = IntegralHelper{Float64}(molecule = M, eri_type=SparseERI)
    rhf = RHF(aoints)
    moints = IntegralHelper(molecule = M, orbitals=rhf.orbitals)
    RMP2(moints, aoints)
end

function RMP2(moints::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,<:AbstractERI,AtomicOrbitals}) where T<:AbstractFloat
    mp_header()
    Fermi.Integrals.compute!(moints, aoints, "OVOV")
    Fermi.Integrals.compute!(moints, aoints, "F")
    RMP2(moints)
end

function RMP2(moints::IntegralHelper{T,E,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,E,AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI}
    mp_header()
    Fermi.Integrals.compute!(moints, aoints, "BOV")
    Fermi.Integrals.compute!(moints, aoints, "F")
    RMP2(moints)
end

function RMP2(ints::IntegralHelper{<:AbstractFloat,<:AbstractERI,<:AbstractRestrictedOrbitals})

    # Check frozen core and inactive virtual
    ndocc = ints.molecule.Nα
    nvir = size(ints.orbitals.C,1) - ndocc
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")

    if core ≥ ndocc
        throw(FermiException("invalid number of frozen orbitals ($core) for $ndocc doubly occupied orbitals"))
    end
    if inac ≥ nvir
        throw(FermiException("invalid number of inactive virtual orbitals ($inac) for $nvir total virtual orbitals"))
    end

    output("  Starting MP2 computation")
    output(" Number of frozen orbitals:             {:d}" , core)
    output(" Number of inactive orbitals:           {:d}" , inac)
    output(" Number of correlated electron pairs:   {:d}", ndocc-core)
    output(" Number of correlated virtual orbitals: {:d}", nvir-inac)
    output(" ⇒ Total number of MP2 amplitudes:      {:d}", (ndocc-core)^2*(nvir-inac)^2)

    output(repeat("-",80))

    # Compute MP2 energy
    Emp2 = RMP2_energy(ints)
    Eref = ints.orbitals.sd_energy

    output("   @Final RMP2 Correlation Energy {:>20.12f} Eₕ", Emp2)
    output("   Reference Energy               {:>20.12f} Eₕ", Eref)
    output("   @Final RMP2 Total Energy       {:>20.12f} Eₕ", Emp2+Eref)
    output(repeat("-",80))

    RMP2(Emp2, Emp2+Eref)
end

function RMP2_energy(ints::IntegralHelper{T, <:AbstractDFERI, RHFOrbitals}) where T<:AbstractFloat
    Bvo = permutedims(ints["BOV"], (1,3,2))
    ϵo = ints["Fii"]
    ϵv = ints["Faa"]

    output(" Computing DF-MP2 Energy!")
    v_size = length(ϵv)
    o_size = length(ϵo)

    # Disable BLAS threading: parallelism happens over occupied index `i` below,
    # each task running single-threaded BLAS (mul!) to avoid oversubscription.
    nt  = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    TWO = T(2.0)
    ONE = one(T)

    # Streaming worker-pool: each task owns its own local Bab buffer and energy
    # accumulator, pulling `i` values off a shared queue, instead of indexing
    # per-thread arrays by Threads.threadid() (unsafe under task migration).
    ntasks = Threads.nthreads()
    requests = Channel{Int}(o_size)
    for i in 1:o_size
        put!(requests, i)
    end
    close(requests)

    results = Channel{T}(ntasks)

    t = @elapsed begin
        @sync for _ in 1:ntasks
            Threads.@spawn begin
            Bab = zeros(T, v_size, v_size)
            ΔMP2_local = zero(T)

            for i in requests
                @views Bi = Bvo[:,:,i]

                for j in i:o_size

                    @views Bj = Bvo[:,:,j]
                    mul!(Bab, transpose(Bi), Bj)

                    eij = ϵo[i] + ϵo[j]
                    E = zero(T)
                    @inbounds for a = eachindex(ϵv)
                        eija = eij - ϵv[a]
                        for b = eachindex(ϵv)
                            D = eija - ϵv[b]
                            E += Bab[a,b] * (TWO * Bab[a,b] - Bab[b,a]) / D
                        end
                    end
                    fac = i !== j ? TWO : ONE
                    ΔMP2_local += fac * E
                end
            end
            put!(results, ΔMP2_local)
            end # spawn
        end # sync

        Emp2 = zero(T)
        for _ in 1:ntasks
            Emp2 += take!(results)
        end
    end # time
    output("Done in {:5.5f} seconds.", t)

    BLAS.set_num_threads(nt)
    return Emp2
end

function RMP2_energy(ints::IntegralHelper{T, Chonky, RHFOrbitals}) where T<:AbstractFloat
    output(" Computing MP2 Energy... ", ending="")
    ovov = ints["OVOV"]
    ϵo = ints["Fii"]
    ϵv = ints["Faa"]

    TWO = T(2.0)

    # Streaming worker-pool: each task keeps a local accumulator and pulls `b`
    # values off a shared queue, instead of indexing a per-thread array by
    # Threads.threadid() (unsafe under task migration / dynamic scheduling).
    ntasks = Threads.nthreads()
    requests = Channel{Int}(length(ϵv))
    for b in eachindex(ϵv)
        put!(requests, b)
    end
    close(requests)

    results = Channel{T}(ntasks)

    t = @elapsed begin
        @sync for _ in 1:ntasks
            Threads.@spawn begin
            ΔMP2_local = zero(T)
            for b in requests
                for a in eachindex(ϵv)
                    @inbounds ϵ_ab = ϵv[a] + ϵv[b]
                    for j in eachindex(ϵo)
                        @inbounds ϵ_abj = ϵo[j] - ϵ_ab
                        for i in eachindex(ϵo)
                            @fastmath @inbounds ΔMP2_local += ovov[i,a,j,b]*(TWO*ovov[i,a,j,b] - ovov[i,b,j,a]) / (ϵo[i]+ϵ_abj)
                        end
                    end
                end
            end
            put!(results, ΔMP2_local)
            end
        end

        Emp2 = zero(T)
        for _ in 1:ntasks
            Emp2 += take!(results)
        end
    end
    output("Done in {:5.5f} s\n", t)
    return Emp2
end

function RMP2_energy(ints::IntegralHelper{T, Chonky, <:AbstractRestrictedOrbitals}) where T<:AbstractFloat
    output(" Computing Iterative MP2 Energy")

    foo = ints["Fij"]
    fvv = ints["Fab"]
    ϵo = ints["Fii"]
    ϵv = ints["Faa"]
    Vovov = ints["OVOV"]
    newT2 = permutedims(Vovov, (1,3,2,4))
    invD = [1.0/(ϵo[i]+ϵo[j]-ϵv[a]-ϵv[b]) for i=eachindex(ϵo), j=eachindex(ϵo), a=eachindex(ϵv), b=eachindex(ϵv)]
    T2 = similar(newT2)

    # Iteration parameters
    ite = 1
    Emp2 = 0.0
    oldE = 0.0
    dE = 1.0
    rms = 1.0
    main_time = 0.0

    e_conv = Options.get("cc_e_conv")
    max_iter = Options.get("cc_max_iter")
    max_rms = Options.get("cc_max_rms")

    output("{:10s}    {: 15s}    {: 12s}    {:12s}    {:10s}","Iteration","CC Energy","ΔE","Max RMS","time (s)")
    while (abs(dE) > e_conv || rms > max_rms)
        if ite > max_iter
            output("\n⚠️  MP2 Equations did not converge in {:1.0d} iterations.", max_iter)
            break
        end
        t = @elapsed begin

            T2 .= newT2
            oldE = Emp2
            TWO = T(2)
            @tensor begin
                newT2[i,j,a,b]  = fvv[c,a]*T2[i,j,c,b]
                newT2[i,j,a,b] -= foo[i,k]*T2[k,j,a,b]
            end
            newT2 .+= permutedims(newT2, (2,1,4,3))
            newT2 .+= permutedims(Vovov, (1,3,2,4))
            newT2 .*= invD

            @tensor begin
                #Energy
                B[l,c,k,d] := TWO*newT2[k,l,c,d]
                B[l,c,k,d] -= newT2[l,k,c,d]
                Emp2 = B[l,c,k,d]*Vovov[k,c,l,d]
            end

            # Compute residues
            rms = sqrt(sum((newT2 .- T2).^2)/length(T2))
        end
        dE = Emp2 - oldE
        main_time += t
        output("    {:<5.0d}    {:<15.10f}    {:>+12.10f}    {:<12.10f}    {:<10.5f}", ite, Emp2, dE, rms, t)
        ite += 1
    end
    output("\nMain MP2 iterations done in {:5.5f} s", main_time)
    output("Average time per iteration {:5.5f}", main_time/(ite-1))

    # Converged?
    conv = false
    if abs(dE) < e_conv && rms < max_rms
        output("\n 🍾 Equations Converged!")
        conv = true
    end

    return Emp2
end
