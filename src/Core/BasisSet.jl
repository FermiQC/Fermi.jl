module GaussianBasis

using Fermi
using Fermi.Error
using Fermi.Geometry
using Fermi.Options
using Fermi.Libcint
using Memoize

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

    output("\n  ⇒ Preparing new Basis Set: {:s}", basis_name)
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
    empty!(memoize_cache(read_basisset))
    output("")
    return BasisSet(mol, basis_name, shells, natm, nbas, nshells, lc_atm, lc_bas, env)
end

end #module