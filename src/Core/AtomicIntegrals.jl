function compute!(I::IntegralHelper, entry::String)

    if entry == "S" 
        compute_S!(I)
    elseif entry == "T" 
        compute_T!(I)
    elseif entry == "V" 
        compute_V!(I)
    elseif entry == "ERI" 
        compute_ERI!(I)
    else
        throw(FermiException("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute_S!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["S"] = FermiMDArray(ao_1e(bs, "overlap", T))
end

function compute_T!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["T"] = FermiMDArray(ao_1e(bs, "kinetic", T))
end

function compute_V!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    bs = I.orbitals.basisset
    I.cache["V"] = FermiMDArray(ao_1e(bs, "nuclear", T))
end

function compute_ERI!(I::IntegralHelper{T, E, AtomicOrbitals}) where {T<:AbstractFloat, E<:AbstractDFERI}

    bs = I.orbitals.basisset
    auxbs = I.eri_type.basisset
    J = FermiMDArray(ao_2e2c(auxbs, T))
    Pqp = FermiMDArray(ao_2e3c(bs, auxbs, T))
    Jh = Array(real(J^(-1/2)))
    @tensor b[Q,p,q] := Pqp[p,q,P]*Jh[P,Q]
    I.cache["ERI"] = b
end

function compute_ERI!(I::IntegralHelper{T, Chonky, AtomicOrbitals}) where T<:AbstractFloat
    bs = I.orbitals.basisset
    I.cache["ERI"] = FermiMDArray(ao_2e4c(bs, T))
end

function compute_ERI!(I::IntegralHelper{T, UniqueERI, AtomicOrbitals}) where T<:AbstractFloat
    bs = I.orbitals.basisset
    I.cache["ERI"] = Fermi.Fermi4SymArray(unique_ao_2e4c(bs, T))
end

function compute_ERI!(I::IntegralHelper{T, SparseERI, AtomicOrbitals}) where T<:AbstractFloat
    bs = I.orbitals.basisset
    I.cache["ERI"] = Fermi.FermiSparse(sparse_ao_2e4c(bs, T)...)
end

function ao_1e(BS::BasisSet, compute::String, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_kin_sph!
    elseif compute == "nuclear"
        libcint_1e! =  cint1e_nuc_sph!
    end

    # Save a list containing the number of primitives for each shell
    num_prim = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]

    # Get slice corresponding to the address in S where the compute chunk goes
    ranges = UnitRange{Int64}[]
    iaccum = 1
    for i = 1:BS.nshells
        push!(ranges, iaccum:(iaccum+ num_prim[i] -1))
        iaccum += num_prim[i]
    end

    # Allocate output array
    out = zeros(T, BS.nbas, BS.nbas)
    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                # Get number of primitive functions
                i_num_prim = num_prim[i]

                # Get slice 
                ri = ranges[i]

                for j in i:BS.nshells

                    # Get number of primitive functions
                    j_num_prim = num_prim[j]

                    # Get slice 
                    rj = ranges[j]

                    # Array where results are written 
                    buf = zeros(Cdouble, i_num_prim*j_num_prim)

                    # Call libcint
                    libcint_1e!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                    # Save results into out
                    out[ri, rj] .= reshape(buf, (i_num_prim, j_num_prim))
                    if j != i
                        out[rj, ri] .= transpose(out[ri, rj])
                    end    
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end

function index2(i::Signed, j::Signed)::Signed
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

# This function is not useful currently, but it will be for direct computations
function index4(i::Signed , j::Signed, k::Signed, l::Signed)::Signed
    return index2(index2(i,j), index2(k,l))
end

# Produces all unique indices ijkl for the two-electron integral
# Note, these indexes start from zero, because we use it in a ccall
function find_indices(nbf::Signed)

    out = NTuple{4,Int16}[]
    N = Int16(nbf - 1)
    ZERO = zero(Int16)

    for i = ZERO:N
        for j = i:N
            for k = ZERO:N
                for l = k:N
                    if index2(i,j) < index2(k,l)
                        continue
                    end
                    push!(out, (i,j,k,l))
                end
            end
        end
    end

    return out
end

function ordered_find_indices(nbf::Signed)

    N = Int16(nbf - 1)
    ZERO = zero(Int16)

    num_ij = Int((nbf^2 - nbf)/2) + nbf
    ij_vals = Array{NTuple{2,Int16}}(undef, num_ij)

    for i = ZERO:N
        for j = i:N
            idx = index2(i,j) + 1
            ij_vals[idx] = (i,j)
        end
    end

    Nelem = Int((num_ij^2 - num_ij)/2) + num_ij 

    out = Array{NTuple{4,Int16}}(undef, Nelem)
    for ij in 1:num_ij
        for kl in 1:ij
            i,j = ij_vals[ij]
            k,l = ij_vals[kl]
            idx = index2(ij-1, kl-1) + 1
            out[idx] = i,j,k,l
        end
    end

    return out
end

function sparse_ao_2e4c(BS::BasisSet, T::DataType = Float64)

    cutoff = Options.get("eri_cutoff")

    # Number of unique integral elements
    N = Int((BS.nbas^2 - BS.nbas)/2) + BS.nbas
    N = Int((N^2 - N)/2) + N

    # Pre allocate output
    out = zeros(T, N)
    indexes = Array{NTuple{4,Int16}}(undef, N)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) - 1 for i = 1:BS.nshells]

    # Unique shell pairs with i < j
    num_ij = Int((BS.nshells^2 - BS.nshells)/2) + BS.nshells

    # Pre allocate array to save ij pairs
    ij_vals = Array{NTuple{2,Int32}}(undef, num_ij)

    # Pre allocate array to save σij, that is the screening parameter for Schwarz 
    σvals = zeros(Float64, num_ij)

    ### Loop thorugh i and j such that i ≤ j. Save each pair into ij_vals and compute √σij for integral screening
    lim = Int32(BS.nshells - 1)
    for i = UnitRange{Int32}(zero(Int32),lim)
        @inbounds begin
        Li2 = lvals[i+1]^2
            for j = UnitRange{Int32}(i, lim)
                Lj2 = lvals[j+1]^2
                buf = zeros(Cdouble, Li2*Lj2)
                idx = index2(i,j) + 1
                ij_vals[idx] = (i,j)
                cint2e_sph!(buf, [i,i,j,j], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                σvals[idx] = √maximum(buf)
            end
        end
    end

    buf_arrays = [zeros(Cdouble, Lmax^4) for _ = 1:Threads.nthreads()]
    
    # i,j,k,l => Shell indexes starting at zero
    # I, J, K, L => AO indexes starting at one
    @sync for ij in eachindex(ij_vals)
        Threads.@spawn begin
        @inbounds begin
            buf = buf_arrays[Threads.threadid()]
            i,j = ij_vals[ij]
            Li, Lj = lvals[i+1], lvals[j+1]
            Lij = Li*Lj
            ioff = ao_offset[i+1]
            joff = ao_offset[j+1]
            for kl in ij:num_ij
                σ = σvals[ij]*σvals[kl]
                if σ < cutoff
                    continue
                end
                k,l = ij_vals[kl]
                Lk, Ll = lvals[k+1], lvals[l+1]
                Lijk = Lij*Lk
                koff = ao_offset[k+1]
                loff = ao_offset[l+1]

                # Compute ERI
                cint2e_sph!(buf, [i,j,k,l], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                ### This block aims to retrieve unique elements within buf and map them to AO indexes
                # is, js, ks, ls are indexes within the shell e.g. for a p shell is = (1, 2, 3)
                # bl, bkl, bjkl are used to map the (i,j,k,l) index into a one-dimensional index for buf
                # That is, get the correct integrals for the AO quartet.
                for ls = 1:Ll
                    L = loff + ls
                    bl = Lijk*(ls-1)
                    for ks = 1:Lk
                        K = koff + ks
                        L < K ? break : nothing

                        # L ≥ K
                        KL = Int(L * (L + 1) / 2) + K                            
                        bkl = Lij*(ks-1) + bl
                        for js = 1:Lj
                            J = joff + js
                            bjkl = Li*(js-1) + bkl
                            for is = 1:Li
                                I = ioff + is
                                J < I ? break : nothing

                                IJ = Int(J * (J + 1) / 2) + I

                                #KL < IJ ? continue : nothing # This restriction does not work... idk why 

                                idx = index2(IJ,KL) + 1
                                out[idx] = buf[is + bjkl]
                                indexes[idx] = (I, J, K, L)
                            end
                        end
                    end
                end
            end
        end #inbounds
        end #spawn
    end #sync
    mask = abs.(out) .> cutoff
    return indexes[mask], out[mask]
end

function unique_ao_2e4c(BS::BasisSet, T::DataType = Float64)

    # Number of unique integral elements
    N = Int((BS.nbas^2 - BS.nbas)/2) + BS.nbas
    N = Int((N^2 - N)/2) + N

    # Pre allocate output
    out = zeros(T, N)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) - 1 for i = 1:BS.nshells]

    # Unique shell pairs with i < j
    num_ij = Int((BS.nshells^2 - BS.nshells)/2) + BS.nshells

    # Pre allocate array to save ij pairs
    ij_vals = Array{NTuple{2,Int32}}(undef, num_ij)

    # Pre allocate array to save σij, that is the screening parameter for Schwarz 
    σvals = zeros(Float64, num_ij)

    ### Loop thorugh i and j such that i ≤ j. Save each pair into ij_vals and compute √σij for integral screening
    lim = Int32(BS.nshells - 1)
    for i = UnitRange{Int32}(zero(Int32),lim)
        @inbounds begin
        Li2 = lvals[i+1]^2
            for j = UnitRange{Int32}(i, lim)
                Lj2 = lvals[j+1]^2
                buf = zeros(Cdouble, Li2*Lj2)
                idx = index2(i,j) + 1
                ij_vals[idx] = (i,j)
                cint2e_sph!(buf, [i,i,j,j], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
                σvals[idx] = √maximum(buf)
            end
        end
    end

    buf_arrays = [zeros(Cdouble, Lmax^4) for _ = 1:Threads.nthreads()]
    
    # i,j,k,l => Shell indexes starting at zero
    # I, J, K, L => AO indexes starting at one
    @sync for ij in eachindex(ij_vals)
        Threads.@spawn begin
        @inbounds begin
            buf = buf_arrays[Threads.threadid()]
            i,j = ij_vals[ij]
            Li, Lj = lvals[i+1], lvals[j+1]
            Lij = Li*Lj
            ioff = ao_offset[i+1]
            joff = ao_offset[j+1]
            for kl in ij:num_ij
                σ = σvals[ij]*σvals[kl]
                if σ < 1e-12
                    continue
                end
                k,l = ij_vals[kl]
                Lk, Ll = lvals[k+1], lvals[l+1]
                Lijk = Lij*Lk
                koff = ao_offset[k+1]
                loff = ao_offset[l+1]

                # Compute ERI
                cint2e_sph!(buf, [i,j,k,l], BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                ### This block aims to retrieve unique elements within buf and map them to AO indexes
                # is, js, ks, ls are indexes within the shell e.g. for a p shell is = (1, 2, 3)
                # bl, bkl, bjkl are used to map the (i,j,k,l) index into a one-dimensional index for buf
                # That is, get the correct integrals for the AO quartet.
                for ls = 1:Ll
                    L = loff + ls
                    bl = Lijk*(ls-1)
                    for ks = 1:Lk
                        K = koff + ks
                        L < K ? break : nothing

                        # L ≥ K
                        KL = Int(L * (L + 1) / 2) + K                            
                        bkl = Lij*(ks-1) + bl
                        for js = 1:Lj
                            J = joff + js
                            bjkl = Li*(js-1) + bkl
                            for is = 1:Li
                                I = ioff + is
                                J < I ? break : nothing

                                IJ = Int(J * (J + 1) / 2) + I

                                #KL < IJ ? continue : nothing # This restriction does not work... idk why 

                                idx = index2(IJ,KL) + 1
                                out[idx] = buf[is + bjkl]
                            end
                        end
                    end
                end
            end
        end #inbounds
        end #spawn
    end #sync
    return out
end

function ao_2e4c(BS::BasisSet, T::DataType = Float64)

    # Save a list containing the number of primitives for each shell
    num_prim = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]

    # Get slice corresponding to the address in S where the compute chunk goes
    ranges = UnitRange{Int64}[]
    iaccum = 1
    for i = 1:BS.nshells
        push!(ranges, iaccum:(iaccum+ num_prim[i] -1))
        iaccum += num_prim[i]
    end

    # Allocate output array
    out = zeros(T, BS.nbas, BS.nbas, BS.nbas, BS.nbas)
    unique_idx = find_indices(BS.nshells)

    @sync for (i,j,k,l) in unique_idx
        Threads.@spawn begin
            @inbounds begin
                # Shift indexes (C starts with 0, Julia 1)
                id, jd, kd, ld = i+1, j+1, k+1, l+1

                # Initialize array for results
                buf = ones(Cdouble, num_prim[id]*num_prim[jd]*num_prim[kd]*num_prim[ld])

                # Compute ERI
                cint2e_sph!(buf, Cint.([i,j,k,l]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                # Move results to output array
                ri, rj, rk, rl = ranges[id], ranges[jd], ranges[kd], ranges[ld]
                out[ri, rj, rk, rl] .= reshape(buf, (num_prim[id], num_prim[jd], num_prim[kd], num_prim[ld]))

                if i != j && k != l && index2(i,j) != index2(k,l)
                    # i,j permutation
                    out[rj, ri, rk, rl] .= permutedims(out[ri, rj, rk, rl], (2,1,3,4))
                    # k,l permutation
                    out[ri, rj, rl, rk] .= permutedims(out[ri, rj, rk, rl], (1,2,4,3))

                    # i,j + k,l permutatiom
                    out[rj, ri, rl, rk] .= permutedims(out[ri, rj, rk, rl], (2,1,4,3))

                    # ij, kl permutation
                    out[rk, rl, ri, rj] .= permutedims(out[ri, rj, rk, rl], (3,4,1,2))
                    # ij, kl + k,l permutation
                    out[rl, rk, ri, rj] .= permutedims(out[ri, rj, rk, rl], (4,3,1,2))
                    # ij, kl + i,j permutation
                    out[rk, rl, rj, ri] .= permutedims(out[ri, rj, rk, rl], (3,4,2,1))
                    # ij, kl + i,j + k,l permutation
                    out[rl, rk, rj, ri] .= permutedims(out[ri, rj, rk, rl], (4,3,2,1))

                elseif k != l && index2(i,j) != index2(k,l)
                    # k,l permutation
                    out[ri, rj, rl, rk] .= permutedims(out[ri, rj, rk, rl], (1,2,4,3))
                    # ij, kl permutation
                    out[rk, rl, ri, rj] .= permutedims(out[ri, rj, rk, rl], (3,4,1,2))
                    # ij, kl + k,l permutation
                    out[rl, rk, ri, rj] .= permutedims(out[ri, rj, rk, rl], (4,3,1,2))

                elseif i != j && index2(i,j) != index2(k,l)
                    # i,j permutation
                    out[rj, ri, rk, rl] .= permutedims(out[ri, rj, rk, rl], (2,1,3,4))

                    # ij, kl permutation
                    out[rk, rl, ri, rj] .= permutedims(out[ri, rj, rk, rl], (3,4,1,2))
                    # ij, kl + i,j permutation
                    out[rk, rl, rj, ri] .= permutedims(out[ri, rj, rk, rl], (3,4,2,1))
        
                elseif i != j && k != l 
                    # i,j permutation
                    out[rj, ri, rk, rl] .= permutedims(out[ri, rj, rk, rl], (2,1,3,4))
                    # k,l permutation
                    out[ri, rj, rl, rk] .= permutedims(out[ri, rj, rk, rl], (1,2,4,3))

                    # i,j + k,l permutatiom
                    out[rj, ri, rl, rk] .= permutedims(out[ri, rj, rk, rl], (2,1,4,3))
                elseif index2(i,j) != index2(k,l) 
                    # ij, kl permutation
                    out[rk, rl, ri, rj] .= permutedims(out[ri, rj, rk, rl], (3,4,1,2))
                end
            end #inbounds
        end #spawn
    end #sync

    return out
end

function ao_2e2c(BS::BasisSet, T::DataType = Float64)

    # Save a list containing the number of primitives for each shell
    num_prim = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]

    # Get slice corresponding to the address in S where the compute chunk goes
    ranges = UnitRange{Int64}[]
    iaccum = 1
    for i = 1:BS.nshells
        push!(ranges, iaccum:(iaccum+ num_prim[i] -1))
        iaccum += num_prim[i]
    end

    # Allocate output array
    out = zeros(T, BS.nbas, BS.nbas)
    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                # Get number of primitive functions
                i_num_prim = num_prim[i]

                # Get slice 
                ri = ranges[i]

                for j in i:BS.nshells

                    # Get number of primitive functions
                    j_num_prim = num_prim[j]

                    # Get slice 
                    rj = ranges[j]

                    # Array where results are written 
                    buf = zeros(Cdouble, i_num_prim*j_num_prim)

                    # Call libcint
                    cint2c2e_sph!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                    # Save results into out
                    out[ri, rj] .= reshape(buf, (i_num_prim, j_num_prim))
                    if j != i
                        out[rj, ri] .= transpose(out[ri, rj])
                    end    
                end
            end #inbounds
        end #spwan
    end #sync
    return out
end

function ao_2e3c(BS::BasisSet, auxBS::BasisSet, T::DataType = Float64)

    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = BS.natoms
    nbas = BS.nbas + auxBS.nbas
    nshells = BS.nshells + auxBS.nshells

    lc_atm = zeros(Cint, natm*ATM_SLOTS)
    lc_bas = zeros(Cint, nshells*BAS_SLOTS)
    env = zeros(Cdouble, length(BS.lc_env) + length(auxBS.lc_env) - 3*natm)

    # Prepare the lc_atom input 
    off = 0
    ib = 0 
    for i = eachindex(BS.molecule.atoms)
        A = BS.molecule.atoms[i]
        # lc_atom has ATM_SLOTS (6) "spaces" for each atom
        # The first one (Z_INDEX) is the atomic number
        lc_atm[1 + ATM_SLOTS*(i-1)] = Cint(A.Z)
        # The second one is the env index address for xyz
        lc_atm[2 + ATM_SLOTS*(i-1)] = off
        env[off+1:off+3] .= A.xyz ./ Fermi.PhysicalConstants.bohr_to_angstrom
        off += 3
        # The remaining 4 slots are zero.
    end

    for i = eachindex(BS.molecule.atoms)
        A = BS.molecule.atoms[i]
        # Prepare the lc_bas input
        for j = eachindex(BS.basis[A])
            B = BS.basis[A][j] 
            Ne = length(B.exp)
            Nc = length(B.coef)
            # lc_bas has BAS_SLOTS for each basis set
            # The first one is the index of the atom starting from 0
            lc_bas[1 + BAS_SLOTS*ib] = i-1
            # The second one is the angular momentum
            lc_bas[2 + BAS_SLOTS*ib] = B.l
            # The third is the number of primitive functions
            lc_bas[3 + BAS_SLOTS*ib] = Nc
            # The fourth is the number of contracted functions
            lc_bas[4 + BAS_SLOTS*ib] = 1
            # The fifth is a κ parameter
            lc_bas[5 + BAS_SLOTS*ib] = 0
            # Sixth is the env index address for exponents
            lc_bas[6 + BAS_SLOTS*ib] = off
            env[off+1:off+Ne] .= B.exp
            off += Ne
            # Seventh is the env index address for contraction coeff
            lc_bas[7 + BAS_SLOTS*ib] = off
            env[off+1:off+Nc] .= B.coef
            off += Nc
            # Eigth, nothing
            ib += 1
        end
    end

    for i = eachindex(auxBS.molecule.atoms)
        A = auxBS.molecule.atoms[i]
        # Prepare the lc_bas input
        for j = eachindex(auxBS.basis[A])
            B = auxBS.basis[A][j] 
            Ne = length(B.exp)
            Nc = length(B.coef)
            # lc_bas has BAS_SLOTS for each basis set
            # The first one is the index of the atom starting from 0
            lc_bas[1 + BAS_SLOTS*ib] = i-1
            # The second one is the angular momentum
            lc_bas[2 + BAS_SLOTS*ib] = B.l
            # The third is the number of primitive functions
            lc_bas[3 + BAS_SLOTS*ib] = Nc
            # The fourth is the number of contracted functions
            lc_bas[4 + BAS_SLOTS*ib] = 1
            # The fifth is a κ parameter
            lc_bas[5 + BAS_SLOTS*ib] = 0
            # Sixth is the env index address for exponents
            lc_bas[6 + BAS_SLOTS*ib] = off
            env[off+1:off+Ne] .= B.exp
            off += Ne
            # Seventh is the env index address for contraction coeff
            lc_bas[7 + BAS_SLOTS*ib] = off
            env[off+1:off+Nc] .= B.coef
            off += Nc
            # Eigth, nothing
            ib += 1
        end
    end

    # Allocate output array
    out = zeros(T, BS.nbas, BS.nbas, auxBS.nbas)

    # Save a list containing the number of primitives for each shell
    num_prim = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Pnum_prim = [Libcint.CINTcgtos_spheric(P-1, auxBS.lc_bas) for P = 1:auxBS.nshells]

    # Get slice corresponding to the address in S where the compute chunk goes
    ranges = UnitRange{Int64}[]
    iaccum = 1
    for i = eachindex(num_prim)
        push!(ranges, iaccum:(iaccum+ num_prim[i] - 1))
        iaccum += num_prim[i]
    end

    Pranges = UnitRange{Int64}[]
    paccum = 1
    for i = eachindex(Pnum_prim)
        push!(Pranges, paccum:(paccum+ Pnum_prim[i] - 1))
        paccum += Pnum_prim[i]
    end

    # icum accumulate the number of primitives in one dimension (i)
    @sync for P in 1:auxBS.nshells
        Threads.@spawn begin
            @inbounds begin
                # Get number of primitive functions
                P_num_prim = Pnum_prim[P]

                # Get slice corresponding to the address in S where the compute chunk goes
                rP = Pranges[P]
                # y accumulate the number of primitives in one dimension (j)
                for i in 1:BS.nshells

                    # Get number of primitive functions
                    i_num_prim = num_prim[i]

                    # Get slice corresponding to the address in S where the compute chunk goes
                    ri = ranges[i]

                    for j in i:BS.nshells

                        j_num_prim = num_prim[j]

                        buf = zeros(Cdouble, i_num_prim*j_num_prim*P_num_prim)

                        rj = ranges[j]

                        # Call libcint
                        cint3c2e_sph!(buf, Cint.([i-1, j-1, P+BS.nshells-1]), lc_atm, natm, lc_bas, nbas, env)

                        # Save results into out
                        out[ri, rj, rP] .= reshape(buf, (i_num_prim, j_num_prim, P_num_prim))
                        if i != j
                            out[rj, ri, rP] .= permutedims(out[ri, rj, rP], (2,1,3))
                        end
                    end
                end
            end #inbounds
        end #spwan
    end #sync
    return out
end

function ao_1e(BS1::BasisSet, BS2::BasisSet, compute::String, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_kin_sph!
    elseif compute == "nuclear"
        libcint_1e! =  cint1e_nuc_sph!
    end

    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = BS1.natoms
    nbas = BS1.nbas + BS2.nbas
    nshells = BS1.nshells + BS2.nshells

    lc_atm = zeros(Cint, natm*ATM_SLOTS)
    lc_bas = zeros(Cint, nshells*BAS_SLOTS)
    env = zeros(Cdouble, length(BS1.lc_env) + length(BS2.lc_env) - 3*natm)

    # Prepare the lc_atom input 
    off = 0
    ib = 0 
    for i = eachindex(BS1.molecule.atoms)
        A = BS1.molecule.atoms[i]
        # lc_atom has ATM_SLOTS (6) "spaces" for each atom
        # The first one (Z_INDEX) is the atomic number
        lc_atm[1 + ATM_SLOTS*(i-1)] = Cint(A.Z)
        # The second one is the env index address for xyz
        lc_atm[2 + ATM_SLOTS*(i-1)] = off
        env[off+1:off+3] .= A.xyz ./ Fermi.PhysicalConstants.bohr_to_angstrom
        off += 3
        # The remaining 4 slots are zero.
    end

    for i = eachindex(BS1.molecule.atoms)
        A = BS1.molecule.atoms[i]
        # Prepare the lc_bas input
        for j = eachindex(BS1.basis[A])
            B = BS1.basis[A][j] 
            Ne = length(B.exp)
            Nc = length(B.coef)
            # lc_bas has BAS_SLOTS for each basis set
            # The first one is the index of the atom starting from 0
            lc_bas[1 + BAS_SLOTS*ib] = i-1
            # The second one is the angular momentum
            lc_bas[2 + BAS_SLOTS*ib] = B.l
            # The third is the number of primitive functions
            lc_bas[3 + BAS_SLOTS*ib] = Nc
            # The fourth is the number of contracted functions
            lc_bas[4 + BAS_SLOTS*ib] = 1
            # The fifth is a κ parameter
            lc_bas[5 + BAS_SLOTS*ib] = 0
            # Sixth is the env index address for exponents
            lc_bas[6 + BAS_SLOTS*ib] = off
            env[off+1:off+Ne] .= B.exp
            off += Ne
            # Seventh is the env index address for contraction coeff
            lc_bas[7 + BAS_SLOTS*ib] = off
            env[off+1:off+Nc] .= B.coef
            off += Nc
            # Eigth, nothing
            ib += 1
        end
    end

    for i = eachindex(BS2.molecule.atoms)
        A = BS2.molecule.atoms[i]
        # Prepare the lc_bas input
        for j = eachindex(BS2.basis[A])
            B = BS2.basis[A][j] 
            Ne = length(B.exp)
            Nc = length(B.coef)
            # lc_bas has BAS_SLOTS for each basis set
            # The first one is the index of the atom starting from 0
            lc_bas[1 + BAS_SLOTS*ib] = i-1
            # The second one is the angular momentum
            lc_bas[2 + BAS_SLOTS*ib] = B.l
            # The third is the number of primitive functions
            lc_bas[3 + BAS_SLOTS*ib] = Nc
            # The fourth is the number of contracted functions
            lc_bas[4 + BAS_SLOTS*ib] = 1
            # The fifth is a κ parameter
            lc_bas[5 + BAS_SLOTS*ib] = 0
            # Sixth is the env index address for exponents
            lc_bas[6 + BAS_SLOTS*ib] = off
            env[off+1:off+Ne] .= B.exp
            off += Ne
            # Seventh is the env index address for contraction coeff
            lc_bas[7 + BAS_SLOTS*ib] = off
            env[off+1:off+Nc] .= B.coef
            off += Nc
            # Eigth, nothing
            ib += 1
        end
    end

    # Allocate output array
    out = zeros(T, BS1.nbas, BS2.nbas)

    # Save a list containing the number of primitives for each shell
    num_prim1 = [Libcint.CINTcgtos_spheric(i-1, BS1.lc_bas) for i = 1:BS1.nshells]
    num_prim2 = [Libcint.CINTcgtos_spheric(i-1, BS2.lc_bas) for i = 1:BS2.nshells]

    # Get slice corresponding to the address in S where the compute chunk goes
    ranges1 = UnitRange{Int64}[]
    iaccum = 1
    for i = eachindex(num_prim1)
        push!(ranges1, iaccum:(iaccum+ num_prim1[i] - 1))
        iaccum += num_prim1[i]
    end

    ranges2 = UnitRange{Int64}[]
    jaccum = 1
    for i = eachindex(num_prim2)
        push!(ranges2, jaccum:(jaccum+ num_prim2[i] - 1))
        jaccum += num_prim2[i]
    end

    # icum accumulate the number of primitives in one dimension (i)
    @sync for i in 1:BS1.nshells
        Threads.@spawn begin
            @inbounds begin
                # Get number of primitive functions
                i_num_prim = num_prim1[i]

                # Get slice corresponding to the address in S where the compute chunk goes
                ri = ranges1[i]
                # y accumulate the number of primitives in one dimension (j)
                for j in 1:BS2.nshells

                    # Get number of primitive functions
                    j_num_prim = num_prim2[j]

                    # Get slice corresponding to the address in S where the compute chunk goes
                    rj = ranges2[j]

                    buf = zeros(Cdouble, i_num_prim*j_num_prim)

                    # Call libcint
                    libcint_1e!(buf, Cint.([i-1, j+BS1.nshells-1]), lc_atm, natm, lc_bas, nbas, env)

                    # Save results into out
                    out[ri, rj] .= reshape(buf, (i_num_prim, j_num_prim))
                end
            end #inbounds
        end #spwan
    end #sync
    return out
end
