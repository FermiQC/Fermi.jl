
function UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc)
    # Calculate energy
    out = Dα .+ Dβ
    out .*= H
    out .+= Fα .* Dα
    out .+= Fβ .* Dβ
    Ee = 0.5*sum(out)
    return(Ee+Vnuc)
end

function build_fock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, Dα, Dβ, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    # Calculate Fock matrix
    H = ints["T"] + ints["V"]
    Fα .= H
    Fβ .= H
    calcJ!(Jα, Dα, ints)
    calcJ!(Jβ, Dβ, ints)
    calcK!(Kα, Dα, ints)
    calcK!(Kβ, Dβ, ints)
    Fα .+= Jα - Kα + Jβ
    Fβ .+= Jβ - Kβ + Jα
end

function calcJ!(J, D, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    # Calculate Coloumb integrals contracted with D
    ERI = ints["ERI"]
    @tensoropt J[i,j] = ERI[i,j,k,l] * D[l,k]
end

function calcK!(K, D, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    # Calculate Exchange integrals contracted with D
    ERI = ints["ERI"]
    @tensoropt K[i,j] = ERI[i,l,k,j] * D[l,k]
end

function build_fock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, Dα, Dβ, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})
    H = ints["T"] + ints["V"]
    Fα .= H
    Fβ .= H
    calcJ!(Jα, Dα, ints)
    calcJ!(Jβ, Dβ, ints)
    calcK!(Kα, Dα, ints)
    calcK!(Kβ, Dβ, ints)
    Fα .+= Jα - Kα + Jβ
    Fβ .+= Jβ - Kβ + Jα
end

function calcJ!(J, D, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})
    # Calculate Coloumb integrals contracted with D
    b = ints["ERI"]
    @tensoropt J[i,j] = b[Q,i,j] * b[Q,k,l] * D[k,l]
end

function calcK!(K, D, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})
    # Calculate Exchange integrals contracted with D
    b = ints["ERI"]
    @tensoropt K[i,j] = b[Q,i,k] * b[Q,j,l] * D[k,l]
end

function buildD!(D, C, N)
    # Build density matrix
    Co = C[:,1:N]
    @tensoropt D[μ, ν] = Co[μ, i] * Co[ν, i]
end

function build_fock!(Fα, Fβ, Jα, Jβ, Kα, Kβ, Dα, Dβ, ints::IntegralHelper{Float64,<:SparseERI,AtomicOrbitals})
    H = ints["T"] + ints["V"]
    Fα .= H
    Fβ .= H
    calcJK!(Jα, Kα, Dα, ints)
    calcJK!(Jβ, Kβ, Dβ, ints)
    Fα .+= Jα - Kα + Jβ
    Fβ .+= Jβ - Kβ + Jα
end

function calcJK!(J, K, D, ints::IntegralHelper{Float64,<:SparseERI,AtomicOrbitals})
    J .= 0.0
    K .= 0.0
    eri_vals = ints["ERI"].data
    idxs = ints["ERI"].indexes

    nbas = ints.orbitals.basisset.nbas
    packed_size = Int((nbas^2 + nbas)/2)

    # Streaming worker-pool: each task owns its own local Jt/Kt buffers and pulls
    # work off a shared queue, instead of indexing a per-thread array by
    # Threads.threadid() (unsafe under task migration / dynamic scheduling).
    # chunksize=1 used to mean one Vector{Int} (and one Channel put!) per single
    # sparse ERI element -- for eri_vals in the millions that's millions of tiny
    # allocations before any Fock-build work starts. Batching amortizes that away
    # while keeping nchunks well above ntasks for load balancing.
    chunksize = clamp(cld(length(eri_vals), 50*Threads.nthreads()), 1, 1000)
    ntasks = Threads.nthreads()

    chunks = [collect(chunk) for chunk in Iterators.partition(eachindex(eri_vals), chunksize)]
    requests = Channel{Vector{Int}}(length(chunks))
    for c in chunks
        put!(requests, c)
    end
    close(requests)

    results = Channel{Tuple{Vector{Float64},Vector{Float64}}}(ntasks)

    @sync begin
        for _ in 1:ntasks
            Threads.@spawn begin
                Jt = zeros(Float64, packed_size)
                Kt = zeros(Float64, packed_size)
                for chunk in requests
                    for z in chunk
                        i,j,k,l = idxs[z] .+ 1
                        ν = eri_vals[z]
                        @inbounds begin
                            @fastmath begin
                            ij = Fermi.index2(i-1,j-1)
                            ik = Fermi.index2(i-1,k-1) + 1
                            il = Fermi.index2(i-1,l-1) + 1

                            jk = Fermi.index2(j-1,k-1) + 1
                            jl = Fermi.index2(j-1,l-1) + 1

                            kl = Fermi.index2(k-1,l-1)

                            γij = i !== j
                            Xik = i === k ? 2 : 1
                            Xjk = j === k ? 2 : 1

                            γkl = k !== l
                            γab = ij !== kl

                            Xil = i === l ? 2 : 1
                            Xjl = j === l ? 2 : 1
                            if γij && γkl && γab
                                # J
                                Jt[ij+1] += 2*D[k,l]*ν
                                Jt[kl+1] += 2*D[i,j]*ν

                                # K
                                Kt[ik] += Xik*D[j,l]*ν
                                Kt[jk] += Xjk*D[i,l]*ν
                                Kt[il] += Xil*D[j,k]*ν
                                Kt[jl] += Xjl*D[i,k]*ν

                            elseif γkl && γab
                                # J
                                Jt[ij+1] += 2*D[k,l]*ν
                                Jt[kl+1] += 1*D[i,j]*ν

                                # K
                                Kt[ik] += Xik*D[j,l]*ν
                                Kt[il] += Xil*D[j,k]*ν
                            elseif γij && γab
                                # J
                                Jt[ij+1] += 1*D[k,l]*ν
                                Jt[kl+1] += 2*D[i,j]*ν

                                # K
                                Kt[ik] += Xik*D[j,l]*ν
                                Kt[jk] += Xjk*D[i,l]*ν

                            elseif γij && γkl

                                # Only possible if i = k and j = l
                                # and i < j ⇒ i < l

                                # J
                                Jt[ij+1] += 2*D[k,l]*ν

                                # K
                                Kt[ik] += D[j,l]*ν
                                Kt[il] += D[j,k]*ν
                                Kt[jl] += D[i,k]*ν
                            elseif γab
                                # J
                                Jt[ij+1] += 1*D[k,l]*ν
                                Jt[kl+1] += 1*D[i,j]*ν
                                # K
                                Kt[ik] += Xik*D[j,l]*ν
                            else
                                Jt[ij+1] += 1*D[k,l]*ν
                                Kt[ik] += D[j,l]*ν
                            end
                            end
                        end
                    end
                end
                put!(results, (Jt, Kt))
            end
        end
    end

    # Reduce values produced by each worker
    Jpacked = zeros(Float64, packed_size)
    Kpacked = zeros(Float64, packed_size)
    for _ in 1:ntasks
        Jt, Kt = take!(results)
        Jpacked .+= Jt
        Kpacked .+= Kt
    end

    for i::Int16 = 1:nbas
        for j::Int16 = i:nbas
            @inbounds begin
                ij = Fermi.index2(i-1,j-1)
                J[i,j] = Jpacked[ij+1]
                K[i,j] = Kpacked[ij+1]

                J[j,i] = J[i,j]
                K[j,i] = K[i,j]
            end
        end
    end
end

function odadamping(D, Ds, F, Fs)
    dD = D - Ds
    s = tr(Fs * dD)
    c = tr((F - Fs) * (dD))
    if c == 0.0
        # dD is identically zero (e.g. an empty spin channel, Nβ=0 -- that
        # density never has anything to build it from, so it stays exactly
        # zero across every iteration): the quadratic damping model's
        # linear and quadratic coefficients (s, c) both vanish, so every λ
        # leaves D and F unchanged. Take the full step (matches the c<=0
        # boundary case below) instead of evaluating -s/(2*c), which would
        # be a literal 0/0 NaN in the branch condition itself.
        λ = 1.0
    elseif c <= -s/(2*c)
        λ = 1.0
    else
        λ = -s/(2*c)
    end
    Fs .= (1-λ)*Fs + λ*F
    Ds .= (1-λ)*Ds + λ*D
    damp = 1-λ
    F .= Fs
    return damp
end

function UHFSpin(Nα , Nβ, Cα, Cβ, S)
    # Calculate S^2 
    smo=transpose(Cα)*S*Cβ
    sz=0.25*(Nα-Nβ)^2
    smo2=smo.^2
    sxy=0.5*(Nα+Nβ)-sum(smo2[1:Nα,1:Nβ])
    spinsquare=sz+sxy
    return(spinsquare)
end
