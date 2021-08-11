function ao_1e(BS::BasisSet, compute::String, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_kin_sph!
    elseif compute == "nuclear"
        libcint_1e! =  cint1e_nuc_sph!
    end

    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    buf_arrays = [zeros(Cdouble, Lmax^4) for _ = 1:Threads.nthreads()]

    @sync for i in 1:BS.nshells
        Threads.@spawn begin
            @inbounds begin
                Li = lvals[i]
                buf = buf_arrays[Threads.threadid()]
                ioff = ao_offset[i]
                for j in i:BS.nshells
                    Lj = lvals[j]
                    joff = ao_offset[j]

                    # Call libcint
                    libcint_1e!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                    # Loop through shell block and save unique elements
                    for js = 1:Lj
                        J = joff + js
                        for is = 1:Li
                            I = ioff + is
                            J < I ? break : nothing
                            out[I,J] = buf[is + Li*(js-1)]
                            out[J,I] = out[I,J]
                        end
                    end
                end
            end #inbounds
        end #spawn
    end #sync
    return out
end

# Needs optimization
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

function grad_ao_1e(BS::BasisSet, iA, compute::String, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ipovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_ipkin_sph!
    # NUCLEAR DOES NOT WORK
    elseif compute == "nuclear"
        libcint_1e! =  cint1e_ipnuc_sph!
    end

    A = BS.molecule.atoms[iA]

    # Pre allocate output
    out = zeros(T, 3, BS.nbas, BS.nbas)

    # Get shell index (s0) where basis of the desired atom start
    s0 = 0
    for a in BS.molecule.atoms
        if a != A
            s0 += length(BS.basis[a])
        else
            break
        end
    end

    # Shell indexes for basis in the atom A
    Ashells = [(s0 + i -1) for i = 1:length(BS.basis[A])]
    notAshells = Int[]
    for i = 1:BS.nshells
        if !((i-1) in Ashells)
            push!(notAshells, i-1)
        end
    end

    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]
    Lmax = maximum(lvals)
    buf = zeros(Cdouble, 3*Lmax^2)
    for i in Ashells
        Li = lvals[i+1]
        ioff = ao_offset[i+1]
        for j in notAshells
            Lj = lvals[j+1]
            joff = ao_offset[j+1]
            # Call libcint
            libcint_1e!(buf, Cint.([j,i]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)
            I = (ioff+1):(ioff+Li)
            J = (joff+1):(joff+Lj)
            out[:,I,J] .= reshape(buf[1:(3*Li*Lj)], 3,Li,Lj)
            out[:,J,I] .= permutedims(out[:,I,J], (1,3,2))
        end
    end
    return out

    # Pre compute a list of angular momentum numbers (l) for each shell
    #lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    #Lmax = maximum(lvals)

    ## Offset list for each shell, used to map shell index to AO index
    #ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    #buf_arrays = [zeros(Cdouble, Lmax^4) for _ = 1:Threads.nthreads()]

    #@sync for i in 1:BS.nshells
    #    Threads.@spawn begin
    #        @inbounds begin
    #            Li = lvals[i]
    #            buf = buf_arrays[Threads.threadid()]
    #            ioff = ao_offset[i]
    #            for j in i:BS.nshells
    #                Lj = lvals[j]
    #                joff = ao_offset[j]

    #                # Call libcint
    #                libcint_1e!(buf, Cint.([i-1,j-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

    #                # Loop through shell block and save unique elements
    #                for js = 1:Lj
    #                    J = joff + js
    #                    for is = 1:Li
    #                        I = ioff + is
    #                        J < I ? break : nothing
    #                        out[I,J] = buf[is + Li*(js-1)]
    #                        out[J,I] = out[I,J]
    #                    end
    #                end
    #            end
    #        end #inbounds
    #    end #spawn
    #end #sync
    #return out
end