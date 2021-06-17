# Auxiliar functions for RHF computations
"""
    function RHF_core_guess(ints::IntegralHelper)

Compute an initial orbital coefficient matrix using the core guess.
"""
function RHF_core_guess(ints::IntegralHelper)
    output("Using Core Guess")
    S = ints["S"]
    Λ = S^(-1/2)
    F = ints["T"] + ints["V"]
    Ft = Λ*F*Λ'

    # Get orbital energies and transformed coefficients
    _, Ct = diagonalize(Ft, hermitian=true)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ
end

"""
    function RHF_gwh_guess(ints::IntegralHelper)

Compute an initial orbital coefficient matrix using the GWH guess.
"""
function RHF_gwh_guess(ints::IntegralHelper)

    # Form GWH guess
    output("Using GWH Guess")
    molecule = ints.molecule
    S = ints["S"]
    d, U = diagonalize(S, sortby = x->1/abs(x))
    Λ = FermiMDArray(Hermitian(S.data)^(-1/2))
    idxs = [abs(d[i]) > 1E-7 for i = eachindex(d)]

    H = real.(ints["T"] + ints["V"])
    ndocc = molecule.Nα
    nvir = size(S,1) - ndocc
    F = similar(S)

    for i = 1:ndocc+nvir
        F[i,i] = H[i,i]
        for j = i+1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ'*F*Λ

    # Get orbital energies and transformed coefficients
    _, Ct = diagonalize(Ft, hermitian=true)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ         
end

"""
    RHFEnergy(D::FermiMDArray{Float64}, H::FermiMDArray{Float64}, F::FermiMDArray{Float64})

Compute RHF energy given a density matrix `D`, Core Hamiltonian `H` and Fock matrix `F`.
"""
function RHFEnergy(D::FermiMDArray{Float64}, H::FermiMDArray{Float64}, F::FermiMDArray{Float64})
    return sum(D .* (H .+ F))
end

"""
    build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ERI::FermiMDArray{Float64,4})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`.
"""
function build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    ERI = ints["ERI"]
    F .= H
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
end

"""
    build_fock!(F::FermiMDArray{Float64,2}, H::FermiMDArray{Float64,2}, D::FermiMDArray{Float64,2}, b::FermiMDArray{Float64,3})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`
approximated by density fitting.
"""
function build_fock!(F::FermiMDArray{Float64,2}, H::FermiMDArray{Float64,2}, D::FermiMDArray{Float64,2}, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})
    b = ints["ERI"]
    F .= H
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
end

function build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ints::IntegralHelper{Float64,UniqueERI,AtomicOrbitals})
    ERI = ints["ERI"].data
    D = D.data
    F .= H
    nbas = ints.orbitals.basisset.nbas

    Farrays = [zeros(Int((nbas^2 + nbas)/2)) for i = 1:Threads.nthreads()]

    @sync for i = 1:nbas
        Threads.@spawn begin
        Ft = Farrays[Threads.threadid()]
        for j = i:nbas
            @inbounds begin
                ij = Fermi.index2(i-1,j-1) 
                γij = i !== j
                for k = 1:nbas
                    Xik = i === k ? 2 : 1
                    Xjk = j === k ? 2 : 1

                    ik = Fermi.index2(i-1,k-1) + 1
                    jk = Fermi.index2(j-1,k-1) + 1
                    for l = k:nbas
                        kl = Fermi.index2(k-1,l-1) 
                        kl < ij ? continue : nothing

                        ν = ERI[Fermi.index2(ij,kl)+1]
                        abs(ν) < 1e-12 ? continue : nothing

                        γkl = k !== l
                        γab = ij !== kl

                        Xil = i === l ? 2 : 1
                        Xjl = j === l ? 2 : 1

                        il = Fermi.index2(i-1,l-1) + 1
                        jl = Fermi.index2(j-1,l-1) + 1

                        if γij && γkl && γab
                            # J
                            Ft[ij+1] += 4*D[k,l]*ν
                            Ft[kl+1] += 4*D[i,j]*ν

                            # K
                            Ft[ik] -= Xik*D[j,l]*ν
                            Ft[jk] -= Xjk*D[i,l]*ν
                            Ft[il] -= Xil*D[j,k]*ν
                            Ft[jl] -= Xjl*D[i,k]*ν

                        elseif γkl && γab
                            # J
                            Ft[ij+1] += 4*D[k,l]*ν
                            Ft[kl+1] += 2*D[i,j]*ν

                            # K
                            Ft[ik] -= Xik*D[j,l]*ν
                            Ft[il] -= Xil*D[j,k]*ν
                        elseif γij && γab
                            # J
                            Ft[ij+1] += 2*D[k,l]*ν
                            Ft[kl+1] += 4*D[i,j]*ν

                            # K
                            Ft[ik] -= Xik*D[j,l]*ν
                            Ft[jk] -= Xjk*D[i,l]*ν

                        elseif γij && γkl

                            # Only possible if i = k and j = l
                            # and i < j ⇒ i < l

                            # J
                            Ft[ij+1] += 4*D[k,l]*ν

                            # K
                            Ft[ik] -= D[j,l]*ν
                            Ft[il] -= D[j,k]*ν
                            Ft[jl] -= D[i,k]*ν
                        elseif γab
                            # J
                            Ft[ij+1] += 2*D[k,l]*ν
                            Ft[kl+1] += 2*D[i,j]*ν
                            # K
                            Ft[ik] -= Xik*D[j,l]*ν
                        else
                            Ft[ij+1] += 2*D[k,l]*ν
                            Ft[ik] -= D[j,l]*ν
                        end
                    end
                end
            end
        end
    end
    end

    for i::Int16 = 1:nbas
        for j::Int16 = i:nbas
            @inbounds begin
                ij = Fermi.index2(i-1,j-1)
                for k = eachindex(Farrays)
                    F[i,j] += Farrays[k][ij+1]
                end
                F[j,i] = F[i,j]
            end
        end
    end
    F
end

function build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ints::IntegralHelper{Float64,SparseERI,AtomicOrbitals})
    ERI = zip(ints["ERI"].indexes, ints["ERI"].data)
    D = D.data
    F .= H
    nbas = ints.orbitals.basisset.nbas

    eri_vals = ints["ERI"].data
    idxs = ints["ERI"].indexes

    Farrays = [zeros(Int((nbas^2 + nbas)/2)) for i = 1:Threads.nthreads()]

    Threads.@threads for z = eachindex(ints["ERI"].data)
        #i0, j0, k0, l0 = idxs[z]
        i,j,k,l = idxs[z] .+ 1
        ν = eri_vals[z]
        Ft = Farrays[Threads.threadid()]
        @inbounds begin
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
                Ft[ij+1] += 4*D[k,l]*ν
                Ft[kl+1] += 4*D[i,j]*ν

                # K
                Ft[ik] -= Xik*D[j,l]*ν
                Ft[jk] -= Xjk*D[i,l]*ν
                Ft[il] -= Xil*D[j,k]*ν
                Ft[jl] -= Xjl*D[i,k]*ν

            elseif γkl && γab
                # J
                Ft[ij+1] += 4*D[k,l]*ν
                Ft[kl+1] += 2*D[i,j]*ν

                # K
                Ft[ik] -= Xik*D[j,l]*ν
                Ft[il] -= Xil*D[j,k]*ν
            elseif γij && γab
                # J
                Ft[ij+1] += 2*D[k,l]*ν
                Ft[kl+1] += 4*D[i,j]*ν

                # K
                Ft[ik] -= Xik*D[j,l]*ν
                Ft[jk] -= Xjk*D[i,l]*ν

            elseif γij && γkl

                # Only possible if i = k and j = l
                # and i < j ⇒ i < l

                # J
                Ft[ij+1] += 4*D[k,l]*ν

                # K
                Ft[ik] -= D[j,l]*ν
                Ft[il] -= D[j,k]*ν
                Ft[jl] -= D[i,k]*ν
            elseif γab
                # J
                Ft[ij+1] += 2*D[k,l]*ν
                Ft[kl+1] += 2*D[i,j]*ν
                # K
                Ft[ik] -= Xik*D[j,l]*ν
            else
                Ft[ij+1] += 2*D[k,l]*ν
                Ft[ik] -= D[j,l]*ν
            end
        end
    end

    # Reduce values produces by each thread
    for i::Int16 = 1:nbas
        for j::Int16 = i:nbas
            @inbounds begin
                ij = Fermi.index2(i-1,j-1)
                for k = eachindex(Farrays)
                    F[i,j] += Farrays[k][ij+1]
                end
                F[j,i] = F[i,j]
            end
        end
    end
    F
end