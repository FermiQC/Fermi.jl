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
    _, Ct = LinearAlgebra.eigen(Symmetric(Ft), sortby=x->x)

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
    d, U = LinearAlgebra.eigen(Symmetric(S), sortby = x->1/abs(x))
    Λ = Hermitian(S)^(-1/2) |> Array
    idxs = [abs(d[i]) > 1E-7 for i = eachindex(d)]

    H = ints["T"] + ints["V"]
    ndocc = molecule.Nα
    nvir = size(S,1) - ndocc
    F = similar(S)

    for i = 1:(ndocc+nvir)
        F[i,i] = H[i,i]
        for j = i+1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ'*F*Λ

    # Get orbital energies and transformed coefficients
    _, Ct = LinearAlgebra.eigen(Symmetric(Ft), sortby=x->x)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ         
end

"""
    RHFEnergy(D::Array{Float64}, H::Array{Float64}, F::Array{Float64})

Compute RHF energy given a density matrix `D`, Core Hamiltonian `H` and Fock matrix `F`.
"""
function RHFEnergy(D::Array{Float64}, H::Array{Float64}, F::Array{Float64})
    return sum(D .* (H .+ F))
end

"""
    build_fock!(F::Array{Float64}, H::Array{Float64}, D::Array{Float64}, ERI::Array{Float64,4})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`.
"""
function build_fock!(F::Array{Float64}, H::Array{Float64}, D::Array{Float64}, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    ERI = ints["ERI"]
    F .= H
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
end

"""
    build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`
approximated by density fitting.
"""
function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})
    b = ints["ERI"]
    F .= H
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
end

function build_fock!(F::Array{Float64}, H::Array{Float64}, D::Array{Float64}, ints::IntegralHelper{Float64,SparseERI,AtomicOrbitals})
    F .= H
    eri_vals = ints["ERI"].data
    idxs = ints["ERI"].indexes

    nbas = ints.orbitals.basisset.nbas

    # --- BEGIN: replaced threaded accumulation with streaming worker-pool ---
    # chunksize=1 used to mean one Vector{Int} (and one Channel put!) per single
    # sparse ERI element -- for eri_vals in the millions that's millions of tiny
    # allocations before any Fock-build work starts. Batching amortizes that away
    # while keeping nchunks well above ntasks for load balancing.
    chunksize = clamp(cld(length(eri_vals), 50*Threads.nthreads()), 1, 1000)
    ntasks = Threads.nthreads()

    # prepare chunked work (vectors of linear indices)
    chunks = [collect(chunk) for chunk in Iterators.partition(eachindex(eri_vals), chunksize)]
    requests = Channel{Vector{Int}}(length(chunks))
    for c in chunks
        put!(requests, c)
    end
    close(requests)

    # bounded channel to receive completed Ft matrices from workers
    results = Channel{Matrix{Float64}}(ntasks)

    @sync begin
        for _ in 1:ntasks
            Threads.@spawn begin
                Ft = zeros(Float64, nbas, nbas)
                for chunk in requests
                    for z in chunk
                        i, j, k, l = idxs[z] .+ 1
                        ν = eri_vals[z]

                        ij = Fermi.index2(i-1, j-1)
                        kl = Fermi.index2(k-1, l-1)

                        γij = i !== j
                        γkl = k !== l
                        γab = ij !== kl

                        Xik = i === k ? 2.0 : 1.0
                        Xjk = j === k ? 2.0 : 1.0
                        Xil = i === l ? 2.0 : 1.0
                        Xjl = j === l ? 2.0 : 1.0

                        if γij && γkl && γab
                            #J
                            Ft[i,j] += 4.0*D[k,l]*ν
                            Ft[k,l] += 4.0*D[i,j]*ν

                            # K
                            Ft[i,k] -= Xik*D[j,l]*ν
                            Ft[j,k] -= Xjk*D[i,l]*ν
                            Ft[i,l] -= Xil*D[j,k]*ν
                            Ft[j,l] -= Xjl*D[i,k]*ν

                        elseif γkl && γab
                            # J
                            Ft[i,j] += 4.0*D[k,l]*ν
                            Ft[k,l] += 2.0*D[i,j]*ν

                            # K
                            Ft[i,k] -= Xik*D[j,l]*ν
                            Ft[i,l] -= Xil*D[j,k]*ν
                        elseif γij && γab
                            # J
                            Ft[i,j] += 2.0*D[k,l]*ν
                            Ft[k,l] += 4.0*D[i,j]*ν

                            # K
                            Ft[i,k] -= Xik*D[j,l]*ν
                            Ft[j,k] -= Xjk*D[i,l]*ν

                        elseif γij && γkl
                            # Only possible if i = k and j = l
                            # and i < j ⇒ i < l

                            # J
                            Ft[i,j] += 4.0*D[k,l]*ν

                            # K
                            Ft[i,k] -= D[j,l]*ν
                            Ft[i,l] -= D[j,k]*ν
                            Ft[j,l] -= D[i,k]*ν
                        elseif γab
                            # J
                            Ft[i,j] += 2.0*D[k,l]*ν
                            Ft[k,l] += 2.0*D[i,j]*ν
                            # K
                            Ft[i,k] -= Xik*D[j,l]*ν
                        else
                            Ft[i,j] += 2.0*D[k,l]*ν
                            Ft[i,k] -= D[j,l]*ν
                        end
                    end
                end
                put!(results, Ft)
            end
        end
    end

    # Consume worker results
    for _ in 1:ntasks
        Ft = take!(results)
        for i in 1:nbas
            F[i,i] += Ft[i,i]
            for j in (i+1):nbas
                F[i,j] += Ft[i,j] + Ft[j,i]
                F[j,i] = F[i,j]
            end
        end
    end
    return F
end