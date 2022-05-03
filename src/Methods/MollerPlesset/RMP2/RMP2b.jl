using LoopVectorization
using LinearAlgebra
using Octavian

function RMP2_canonical_energy(ints::IntegralHelper{T,RIFIT,O}, Alg::RMP2b) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
    Bvo = permutedims(ints["BOV"], (1,3,2))
    ϵo = ints["Fii"]
    ϵv = ints["Faa"]

    output(" Computing DF-MP2 Energy")
    output(" - Contraction engine: Octavian")
    v_size = length(ϵv)
    o_size = length(ϵo)

    # Pre-allocating arrays for threads
    BABs = [zeros(T, v_size, v_size) for _ = 1:Threads.nthreads()]

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
                matmul_serial!(Bab, transpose(Bi), Bj)

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

    return Emp2
end