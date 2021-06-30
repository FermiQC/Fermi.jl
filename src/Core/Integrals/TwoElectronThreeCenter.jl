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
    #N = (BS.nbas^2 + BS.nbas) >> 1
    #out = zeros(T, N, auxBS.nbas)
    out = zeros(T, BS.nbas, BS.nbas, auxBS.nbas)

    # Save a list containing the angular momentum number for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Plvals = [Libcint.CINTcgtos_spheric(P-1, auxBS.lc_bas) for P = 1:auxBS.nshells]
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]


    buf_arrays = [zeros(Cdouble, maximum(lvals)^2*maximum(Plvals)) for _ = 1:Threads.nthreads()]

    @sync for p in 1:auxBS.nshells
        Threads.@spawn begin
            @inbounds begin
                buf = buf_arrays[Threads.threadid()]
                Lp = Plvals[p]
                poff = sum(Plvals[1:(p-1)])
                Bp = p + BS.nshells - 1
                for i in 1:BS.nshells
                    Li = lvals[i]
                    ioff = ao_offset[i]
                    for j in i:BS.nshells
                        Lj = lvals[j]
                        joff = ao_offset[j]

                        # Call libcint
                        cint3c2e_sph!(buf, Cint.([i-1, j-1, Bp]), lc_atm, natm, lc_bas, nbas, env)

                        # Loop through shell block and save unique elements
                        for ps = 1:Lp
                            P = poff + ps
                            for js = 1:Lj
                                J = joff + js
                                for is = 1:Li
                                    I = ioff + is
                                    J < I ? break : nothing
                                    out[I,J,P] = buf[is + Li*(js-1) + Li*Lj*(ps-1)]
                                    out[J,I,P] = out[I,J,P]
                                end
                            end
                        end
                    end
                end
            end #inbounds
        end #spwan
    end #sync
    return out
end