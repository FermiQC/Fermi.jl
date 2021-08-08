
function UHFEnergy(H, Dα, Dβ, Fα, Fβ, Vnuc)
    # Calculate energy
    @tensoropt A[:] := 0.5 * (H[i,j] * (Dα[i,j] + Dβ[i,j]) + Fα[i,j] * Dα[i,j] + Fβ[i,j] * Dβ[i,j])
    Ee = A[1]
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
    D = D.data
    J .= 0.0
    K .= 0.0
    eri_vals = ints["ERI"].data
    idxs = ints["ERI"].indexes

    nbas = ints.orbitals.basisset.nbas
    Jarrays = [zeros(Int((nbas^2 + nbas)/2)) for i = 1:Threads.nthreads()]
    Karrays = [zeros(Int((nbas^2 + nbas)/2)) for i = 1:Threads.nthreads()]

    Threads.@threads for z = eachindex(eri_vals)
        i,j,k,l = idxs[z] .+ 1
        ν = eri_vals[z]
        Jt = Jarrays[Threads.threadid()]
        Kt = Karrays[Threads.threadid()]
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

    # Reduce values produces by each thread
    for i::Int16 = 1:nbas
        for j::Int16 = i:nbas
            @inbounds begin
                ij = Fermi.index2(i-1,j-1)
                for k = eachindex(Jarrays)
                    J[i,j] += Jarrays[k][ij+1]
                end

                for k = eachindex(Karrays)
                    K[i,j] += Karrays[k][ij+1]
                end

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
    if c <= -s/(2*c)
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