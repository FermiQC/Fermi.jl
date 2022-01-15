using LoopVectorization
using LinearAlgebra

function RMP2(Alg::RMP2Algorithm)
    aoints = IntegralHelper{Float64}()
    rhf = RHF(aoints)

    if typeof(aoints.eri_type) === JKFIT || Options.get("precision") == "single"
        aoints = IntegralHelper(eri_type=RIFIT())
    end
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RMP2(moints, aoints, Alg)
end

function RMP2(aoints::IntegralHelper{Float64,E,AtomicOrbitals}, Alg::RMP2Algorithm) where E <: AbstractERI
    rhf = RHF(aoints)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RMP2(moints, aoints, Alg)
end

function RMP2(O::AbstractRestrictedOrbitals, Alg::RMP2Algorithm)
    moints = IntegralHelper(orbitals=O)
    aoints = IntegralHelper(eri_type=moints.eri_type)
    RMP2(moints, aoints, Alg)
end

function RMP2(rhf::RHF, Alg::RMP2Algorithm)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    aoints = IntegralHelper(eri_type=moints.eri_type)
    RMP2(moints, aoints, Alg)
end

function RMP2(M::Molecule, Alg::RMP2Algorithm)
    aoints = IntegralHelper{Float64}(molecule = M, eri_type=SparseERI)
    rhf = RHF(aoints)
    moints = IntegralHelper(molecule = M, orbitals=rhf.orbitals)
    RMP2(moints, aoints, Alg)
end

function RMP2(moints::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,<:AbstractERI,AtomicOrbitals}, Alg::RMP2Algorithm) where T<:AbstractFloat
    mp_header()
    Fermi.Integrals.compute!(moints, aoints, "OVOV")
    Fermi.Integrals.compute!(moints, aoints, "F")
    RMP2(moints, Alg)
end

function RMP2(moints::IntegralHelper{T,E,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,E,AtomicOrbitals}, Alg::RMP2Algorithm) where {T<:AbstractFloat, E<:AbstractDFERI}
    mp_header()
    Fermi.Integrals.compute!(moints, aoints, "BOV")
    Fermi.Integrals.compute!(moints, aoints, "F")
    RMP2(moints, Alg)
end

function RMP2(ints::IntegralHelper{<:AbstractFloat,<:AbstractERI,<:AbstractRestrictedOrbitals}, Alg::RMP2Algorithm)

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
    Emp2 = RMP2_energy(ints, Alg)
    Eref = ints.orbitals.sd_energy

    output("   @Final RMP2 Correlation Energy {:>20.12f} Eₕ", Emp2)
    output("   Reference Energy               {:>20.12f} Eₕ", Eref)
    output("   @Final RMP2 Total Energy       {:>20.12f} Eₕ", Emp2+Eref)
    output(repeat("-",80))

    RMP2(Emp2, Emp2+Eref)
end

function RMP2_energy(ints::IntegralHelper{T,E,O}, Alg::RMP2Algorithm) where {T<:AbstractFloat, E<:AbstractERI,O<:AbstractRestrictedOrbitals}

    Emp2 = RMP2_nonrhf_energy(ints, Alg)
    Emp2 += RMP2_canonical_energy(ints, Alg)

    return Emp2
end

function RMP2_energy(ints::IntegralHelper{<:AbstractFloat, <:AbstractERI, RHFOrbitals}, Alg::RMP2Algorithm) 
    Emp2 = RMP2_canonical_energy(ints, Alg)
    return Emp2
end

function RMP2_nonrhf_energy(ints::IntegralHelper{<:AbstractFloat, <:AbstractERI, <:AbstractRestrictedOrbitals}, Alg::RMP2Algorithm)

    ϵo = ints["Fii"]
    ϵv = ints["Faa"]

    t = @elapsed begin
    output("Computing non-RHF contribution to the MP2 energy... ", ending="")

    Dia = [ϵo[i]-ϵv[a] for i = eachindex(ϵo), a = eachindex(ϵv)]
    Fia = ints["Fia"]
    Emp2 = sum(transpose(Fia)*(Fia ./ Dia))
    end
    output("Done in {:3.5f} seconds.", t)

    return Emp2
end

function RMP2_canonical_energy(ints::IntegralHelper{T, <:AbstractDFERI, <:AbstractRestrictedOrbitals}, Alg::RMP2Algorithm) where T<:AbstractFloat
    Bvo = permutedims(ints["BOV"].data, (1,3,2))
    ϵo = ints["Fii"]
    ϵv = ints["Faa"]

    output(" Computing DF-MP2 Energy!")
    v_size = length(ϵv)
    o_size = length(ϵo)

    # Pre-allocating arrays for threads
    BABs = [zeros(T, v_size, v_size) for _ = 1:Threads.nthreads()]

    # Disable BLAS threading
    nt  = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    # Vector containing the energy contribution computed by each thread
    ΔMP2s = zeros(T,Threads.nthreads())
    TWO = T(2.0)
    ONE = one(T)
    t = @elapsed begin
        @sync for i in 1:o_size
            Threads.@spawn begin
            id = Threads.threadid()

            Bab = BABs[id]
            @views Bi = Bvo[:,:,i]

            for j in i:o_size

                @views Bj = Bvo[:,:,j]
                mul!(Bab, transpose(Bi), Bj)

                eij = ϵo[i] + ϵo[j]
                E = zero(T)
                @turbo for a = eachindex(ϵv)
                    eija = eij - ϵv[a]
                    for b = eachindex(ϵv)
                        D = eija - ϵv[b]
                        E += Bab[a,b] * (TWO * Bab[a,b] - Bab[b,a]) / D
                    end
                end
                fac = i !== j ? TWO : ONE 
                ΔMP2s[id] += fac * E
            end
        end # spawn
        end # sync
        Emp2 = sum(ΔMP2s)
    end # time
    output("Done in {:5.5f} seconds.", t)

    BLAS.set_num_threads(nt)
    return Emp2
end

function RMP2_canonical_energy(ints::IntegralHelper{T, Chonky, <:AbstractRestrictedOrbitals}, Alg::RMP2Algorithm) where T<:AbstractFloat
    output(" Computing MP2 Energy... ", ending="")
    ovov = ints["OVOV"]
    ϵo = ints["Fii"]
    ϵv = ints["Faa"]

    ΔMP2s = zeros(T,Threads.nthreads())
    TWO = T(2.0)
    t = @elapsed begin
        @sync for b in eachindex(ϵv)
            Threads.@spawn begin
            id = Threads.threadid()
            for a in eachindex(ϵv)
                @inbounds ϵ_ab = ϵv[a] + ϵv[b]
                for j in eachindex(ϵo)
                    @inbounds ϵ_abj = ϵo[j] - ϵ_ab
                    for i in eachindex(ϵo)
                        @fastmath @inbounds ΔMP2s[id] += ovov[i,a,j,b]*(TWO*ovov[i,a,j,b] - ovov[i,b,j,a]) / (ϵo[i]+ϵ_abj)
                    end
                end
            end
            end
        end
    end
    output("Done in {:5.5f} s\n", t)
    return sum(ΔMP2s)
end
