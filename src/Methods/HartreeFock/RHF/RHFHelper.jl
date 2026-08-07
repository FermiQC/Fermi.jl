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
    chunksize = 1
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

"""
    shell_to_atom_index(bset, sh)

Map a basis 'shell' object `sh` to the 1-based atom index in the given
`bset` (GaussianBasis.BasisSet). The function accepts several common
representations that shell objects use to indicate their owning atom:

- an integer index (already 1-based or 0-based depending on provider)
- a field named `:atom` or `:atom_index` containing a numeric index
- a field named `:center` containing either an Atom-like object or a
  coordinate vector

"""
function shell_to_atom_index(bset, sh)
    fn = fieldnames(typeof(sh))
    if :atom in fn
        val = getfield(sh, :atom)
    elseif :atom_index in fn
        val = getfield(sh, :atom_index)
    elseif :center in fn
        val = getfield(sh, :center)
    else
        throw(ErrorException("Cannot determine shell->atom mapping. Shell fields: $(fn)"))
    end

    # If value is already an integer index, return it
    if isa(val, Integer)
        return Int(val)
    end

    # If value is an Atom-like object, try to find its index in bset.atoms
    for (i,a) in enumerate(bset.atoms)
        try
            if a === val
                return i
            end
            if hasfield(typeof(a), :Z) && hasfield(typeof(val), :Z)
                if getfield(a, :Z) == getfield(val, :Z) && isapprox(getfield(a, :xyz), getfield(val, :xyz); atol=1e-8)
                    return i
                end
            end
            if hasfield(typeof(a), :AtomicSymbol) && hasfield(typeof(val), :AtomicSymbol)
                if getfield(a, :AtomicSymbol) == getfield(val, :AtomicSymbol) && isapprox(getfield(a, :xyz), getfield(val, :xyz); atol=1e-8)
                    return i
                end
            end
        catch
            # ignore and continue
        end
    end

    # If the value is a coordinate vector, match by coordinates
    if isa(val, AbstractVector)
        for (i,a) in enumerate(bset.atoms)
            if isapprox(getfield(a, :xyz), val; atol=1e-8)
                return i
            end
        end
    end

    throw(ErrorException("Failed to map shell to atom index for shell value: $(val)"))
end

## This is the old version of build_fock! that uses threaded accumulation into thread-specific arrays, 
# which is no longer recommended. Left here a token of admiraton for the times where threading was simpler to implement and understand RIP.
function old_build_fock!(F::Array{Float64}, H::Array{Float64}, D::Array{Float64}, ints::IntegralHelper{Float64,SparseERI,AtomicOrbitals})
    F .= H
    eri_vals = ints["ERI"].data
    idxs = ints["ERI"].indexes

    nbas = ints.orbitals.basisset.nbas
    Farrays = [zeros(nbas, nbas) for i = 1:Threads.maxthreadid()]

    Threads.@threads :static for z = eachindex(eri_vals)
        i,j,k,l = idxs[z] .+ 1
        ν = eri_vals[z]
        Ft = Farrays[Threads.threadid()]
        ij = Fermi.index2(i-1,j-1) 
        kl = Fermi.index2(k-1,l-1) 

        # Logical auxiliar: γpq (whether p and q are different) Xpq = δpq + 1
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

    Fred = sum(Farrays)
    Threads.@threads :static for i = 1:nbas
        F[i,i] += Fred[i,i] 
        for j = (i+1):nbas
            F[i,j] += Fred[i,j] + Fred[j,i]
            F[j,i] = F[i,j]
        end
    end
    F
end