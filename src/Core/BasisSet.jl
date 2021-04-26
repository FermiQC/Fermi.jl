module GaussianBasis

using Fermi.Geometry
using Fermi.Options
using Fermi.Libcint

import Base: collect, +

export BasisSet, BasisFunction

include("BasisParser.jl")

struct BasisFunction
    l::Cint
    coef::Array{Cdouble,1}
    exp::Array{Cdouble,1}
end

struct BasisSet
    molecule::Molecule
    basis_name::String
    basis::Dict{Atom,Array{BasisFunction,1}}
    natoms::Cint
    nbas::Cint
    nshells::Int64
    lc_atoms::Array{Cint,1}
    lc_bas::Array{Cint, 1}
    lc_env::Array{Cdouble,1}
end

function gto_norm(n::Signed, a::AbstractFloat)
   # normalization factor of function rⁿ exp(-ar²)
    s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
    return √s
end

function normalize_basisfunction!(B::BasisFunction)
    for i = eachindex(B.coef)
        B.coef[i] *= gto_norm(B.l, B.exp[i])
    end
end

function BasisSet()
    mol = Molecule()
    bname = Options.get("basis")
    BasisSet(mol, bname)
end

function BasisSet(mol::Molecule, basis_name::String)

    atoms = mol.atoms
    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = length(atoms)
    nbas = 0
    nshells = 0
    nexps = 0
    nprims = 0
    shells = Dict{Atom, Array{BasisFunction,1}}()

    for A in atoms
        basis = read_basisset(basis_name, A.AtomicSymbol)
        for b in basis
            nshells += 1
            nbas += 2*b.l + 1
            nexps += length(b.exp)
            nprims += length(b.coef)
            normalize_basisfunction!(b)
        end
        shells[A] = basis
    end

    lc_atm = zeros(Cint, natm*ATM_SLOTS)
    lc_bas = zeros(Cint, nshells*BAS_SLOTS)
    env = zeros(Cdouble, 3*natm+nexps+nprims)

    # Prepare the lc_atom input 
    off = 0
    ib = 0 
    for i = eachindex(atoms)
        A = atoms[i]
        # lc_atom has ATM_SLOTS (6) "spaces" for each atom
        # The first one (Z_INDEX) is the atomic number
        lc_atm[1 + ATM_SLOTS*(i-1)] = Cint(A.Z)
        # The second one is the env index address for xyz
        lc_atm[2 + ATM_SLOTS*(i-1)] = off
        env[off+1:off+3] .= A.xyz ./ Fermi.PhysicalConstants.bohr_to_angstrom
        off += 3
        # The remaining 4 slots are zero.

        # Prepare the lc_bas input
        for j = eachindex(shells[A])
            B = shells[A][j] 
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
    return BasisSet(mol, basis_name, shells, natm, nbas, nshells, lc_atm, lc_bas, env)
end

function collect(BS::BasisSet)
    out = BasisFunction[]
    for A in BS.molecule.atoms
        push!(out, BS.basis[A]...)
    end
    return out
end

function +(B1::BasisSet, B2::BasisSet)


    atoms = B1.molecule.atoms


    return BasisSet(B1.molecule, basis_name, B1.basis, natm, nbas, nshells, lc_atom, lc_bas, lc_env)
end

function ao_1e(BS::BasisSet, compute::String)

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
    out = zeros(Cdouble, BS.nbas, BS.nbas)
    @sync for i in 1:BS.nshells
    Threads.@spawn begin

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
    end
    end
    return out
end

function ao_2e4c(BS::BasisSet)

    # x accumulate the number of primitives in one dimension (i)
    icum = 1
    # Allocate output array
    out = zeros(Cdouble, BS.nbas, BS.nbas, BS.nbas, BS.nbas)
    for i in 1:BS.nshells

        # Get number of primitive functions
        i_num_prim = Libcint.CINTcgtos_spheric(i-1, BS.lc_bas)

        # Get slice corresponding to the address in S where the compute chunk goes
        ri = icum:(icum+i_num_prim-1)
        icum += i_num_prim
        # y accumulate the number of primitives in one dimension (j)
        jcum = 1
        for j in 1:BS.nshells

            # Get number of primitive functions
            j_num_prim = Libcint.CINTcgtos_spheric(j-1, BS.lc_bas)

            # Get slice corresponding to the address in S where the compute chunk goes
            rj = jcum:(jcum+j_num_prim-1)
            jcum += j_num_prim

            kcum = 1
            for k in 1:BS.nshells
                # Get number of primitive functions
                k_num_prim = Libcint.CINTcgtos_spheric(k-1, BS.lc_bas)

                # Get slice corresponding to the address in S where the compute chunk goes
                rk = kcum:(kcum+k_num_prim-1)
                kcum += k_num_prim

                lcum = 1
                for l in 1:BS.nshells
                    # Get number of primitive functions
                    l_num_prim = Libcint.CINTcgtos_spheric(l-1, BS.lc_bas)

                    # Get slice corresponding to the address in S where the compute chunk goes
                    rl = lcum:(lcum+l_num_prim -1)
                    lcum += l_num_prim
                    # Array where results are written 
                    buf = zeros(Cdouble, i_num_prim*j_num_prim*k_num_prim*l_num_prim)

                    # Call libcint
                    cint2e_sph!(buf, Cint.([i-1,j-1,k-1,l-1]), BS.lc_atoms, BS.natoms, BS.lc_bas, BS.nbas, BS.lc_env)

                    # Save results into out
                    out[ri, rj, rk, rl] .= reshape(buf, (i_num_prim, j_num_prim, k_num_prim, l_num_prim))
                end
            end
        end
    end
    return out
end

function ao_2e2c(BS::BasisSet)

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
    out = zeros(Cdouble, BS.nbas, BS.nbas)
    @sync for i in 1:BS.nshells
    Threads.@spawn begin

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
    end
    end
    return out
end

function ao_2e3c(BS::BasisSet, auxBS::BasisSet)

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
    out = zeros(Cdouble, BS.nbas, BS.nbas, auxBS.nbas)

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
    Pcum = 1
    @sync for P in 1:auxBS.nshells
    Threads.@spawn begin

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
    end
    end
    return out
end

end #module