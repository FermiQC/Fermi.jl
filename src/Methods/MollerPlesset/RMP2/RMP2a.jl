using TensorOperations
using LinearAlgebra

function RMP2(ints::IntegralHelper{T,E,O}, Alg::RMP2a) where {T<:AbstractFloat,E<:AbstractERI,O<:AbstractRestrictedOrbitals}
    mp_header()

    # Check frozen core and inactive virtual
    ndocc = ints.molecule.Nα
    nvir = size(ints.orbitals.C,1) - ndocc
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")

    if core ≥ ndocc
        throw(InvalidFermiOption("invalid number of frozen orbitals ($core) for $ndocc doubly occupied orbitals"))
    end
    if inac ≥ nvir
        throw(InvalidFermiOption("invalid number of inactive virtual orbitals ($inac) for $nvir total virtual orbitals"))
    end

    # Collect RHF information

    output("  Starting MP2 computation")
    output(" Number of frozen orbitals:             {:d}" , core)
    output(" Number of inactive orbitals:           {:d}" , inac)
    output(" Number of correlated electron pairs:   {:d}", ndocc-core)
    output(" Number of correlated virtual orbitals: {:d}", nvir-inac)
    output(" ⇒ Total number of MP2 amplitudes:      {:d}", (ndocc-core)^2*(nvir-inac)^2)

    output(repeat("-",80))

    # Compute MP2 energy
    t = @elapsed Emp2 = RMP2_energy(ints, Alg)
    Eref = Integrals.reference_energy(ints)

    output("   @Final RMP2 Correlation Energy {:>20.12f} Eₕ", Emp2)
    output("   Reference Energy               {:>20.12f} Eₕ", Eref)
    output("   @Final RMP2 Total Energy       {:>20.12f} Eₕ", Emp2+Eref)
    output(repeat("-",80))

    RMP2(Emp2, Emp2+Eref)
end

function RMP2_energy(ints::IntegralHelper{T,E,O}, Alg::RMP2a) where {T<:AbstractFloat, E<:AbstractERI,O<:AbstractRestrictedOrbitals}
    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = ints.molecule.Nα
    nbf = size(ints.orbitals.C,1)
    ϵo = ints["Fd"][(core+1):ndocc]
    ϵv = ints["Fd"][(ndocc+1):(nbf-inac)]

    Emp2 = RMP2_nonrhf_energy(ints, ϵo, ϵv)
    Emp2 += RMP2_rhf_energy(ints, ϵo, ϵv)

    return Emp2
end

function RMP2_energy(ints::IntegralHelper{T,E,RHFOrbitals}, Alg::RMP2a) where {T<:AbstractFloat, E<:AbstractERI}
    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = ints.molecule.Nα
    nbf = size(ints.orbitals.C,1)
    ϵo = ints["Fd"][(core+1):ndocc]
    ϵv = ints["Fd"][(ndocc+1):(nbf-inac)]

    Emp2 = RMP2_rhf_energy(ints, ϵo, ϵv)

    return Emp2
end

function RMP2_nonrhf_energy(ints::IntegralHelper{T,E,O}, ϵo::AbstractArray{T,1}, ϵv::AbstractArray{T,1}) where {T<:AbstractFloat,E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    t = @elapsed begin
    output("Computing non-RHF contribution to the MP2 energy... ", ending="")

    Dia = [ϵo[i]-ϵv[a] for i = eachindex(ϵo), a = eachindex(ϵv)]
    Fia = ints["Fia"]
    Emp2 = sum(transpose(Fia)*(Fia ./ Dia))
    end
    output("Done in {:3.5f} seconds.", t)

    return Emp2
end

function RMP2_rhf_energy(ints::IntegralHelper{T,RIFIT,O}, ϵo::AbstractArray{T,1}, ϵv::AbstractArray{T,1}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
    Bov = ints["BOV"].data

    output(" Computing DF-MP2 Energy... ", ending="")
    v_size = length(ϵv)
    o_size = length(ϵo)

    # Pre-allocating arrays for threads
    BABs = [zeros(T, v_size, v_size) for i=1:Threads.nthreads()]
    BBAs = [zeros(T, v_size, v_size) for i=1:Threads.nthreads()]

    # Vector containing the energy contribution computed by each thread
    ΔMP2s = zeros(T,Threads.nthreads())
    t = @elapsed begin
        @sync for i in 1:o_size
            Threads.@spawn begin
            id = Threads.threadid()
            @views Bi = Bov[:,i,:]
            BAB = BABs[id]
            BBA = BBAs[id]
            for j in i:o_size
                if i != j
                    fac = T(2)
                else
                    fac = one(T)
                end
                @views Bj = Bov[:,j,:]

                @tensor BAB[a,b] = Bi[Q,a]*Bj[Q,b]
                transpose!(BBA,BAB)

                eij = ϵo[i] + ϵo[j]
                dmp2_1 = zero(T)
                dmp2_2 = zero(T) 
                @fastmath for b in 1:v_size
                    @inbounds eij_b = eij - ϵv[b]
                    for a in 1:v_size
                        @inbounds begin 
                            iajb = BAB[a,b]
                            d = iajb/(eij_b - ϵv[a])
                            dmp2_1 += d*iajb
                            dmp2_2 += d*BBA[a,b]
                        end
                    end
                end
                dmp2 = T(2)*dmp2_1 - dmp2_2
                ΔMP2s[id] += fac*dmp2
            end
            end
        end
        Emp2 = sum(ΔMP2s)
    end
    output("Done in {:5.5f} seconds.", t)
    return Emp2
end

function RMP2_rhf_energy(ints::IntegralHelper{T,Chonky,O}, ϵo::AbstractArray{T,1}, ϵv::AbstractArray{T,1}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
    output(" Computing MP2 Energy... ", ending="")
    ovov = ints["OVOV"]

    ΔMP2s = zeros(T,Threads.nthreads())
    t = @elapsed begin
        Emp2 = zero(T)
        @sync for b in eachindex(ϵv)
            Threads.@spawn begin
            id = Threads.threadid()
            for a in eachindex(ϵv)
                ϵ_ab = ϵv[a] + ϵv[b]
                for j in eachindex(ϵo)
                    ϵ_abj = ϵo[j] - ϵ_ab
                    @fastmath for i in eachindex(ϵo)
                        @inbounds ΔMP2s[id] += ovov[i,a,j,b]*(T(2)*ovov[i,a,j,b] - ovov[i,b,j,a]) /
                        (ϵo[i]+ϵ_abj)
                    end
                end
            end
            end
        end
    end
    output("Done in {:5.5f} s\n", t)
    return sum(ΔMP2s)
end
