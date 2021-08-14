"""
    Fermi.GaussianBasis

Module handling basis functions and basis set.
"""
module GaussianBasis

using Fermi
using Fermi.Options
using Memoize
using Formatting
using Molecules

import Base: getindex, show

export BasisSet, BasisFunction

include("BasisParser.jl")

@doc raw"""
    Fermi.GaussianBasis.BasisFunction

Object representing a shell of Gaussian basis functions composed of ``N`` primitives: 

```math
\chi_{l,m} = \sum_i^N C_i r^{l}e^{-\zeta_i r^2} Y_{l,m}(\theta,\phi)
```

# Fields

| Name    | Type | Description |
|:--------|:------|:-----------------------------------------------------------|
|`l`      |`Int32`            | Angular momentum number (e.g. 0, 1, 2 for S, P, D...)      |
|`coef`   |`Array{Float64,1}` | Array with coefficients (``C_i``) for each primitive       |
|`exp`    |`Array{Float64,1}` | Array with exponents (``\zeta_i``) for each primitive      |

# Examples

```julia
julia> bf = Fermi.GaussianBasis.BasisFunction(1, [1/√2, 1/√2], [5.0, 1.2])
P shell with 3 basis built from 2 primitive gaussians

χ₁₋₁ =    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₋₁⋅r¹⋅exp(-1.2⋅r²)

χ₁₀  =    0.7071067812⋅Y₁₀⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₀⋅r¹⋅exp(-1.2⋅r²)

χ₁₁  =    0.7071067812⋅Y₁₁⋅r¹⋅exp(-5.0⋅r²)
     +    0.7071067812⋅Y₁₁⋅r¹⋅exp(-1.2⋅r²)
```
"""
struct BasisFunction
    l::Cint
    coef::Array{Cdouble,1}
    exp::Array{Cdouble,1}
end

"""
    Fermi.GaussianBasis.BasisSet

Object holding a set of BasisFunction objects associated with each atom in a molecule.

# Fields

| Name        | Type                               |   Description |
|:------------|:-----------------------------------|:-----------------------------------------------------------|
|`molecule`   | `Molecule`                         | Fermi.Molecule object  |
|`basis_name` | `String`                           | String holding the basis set name  |
|`basis`      | `Dict{Atom,Array{BasisFunction,1}}`| A dictionary that maps Atom objects to an Array of BasisFunction  |
|`natoms`     | `Int32`                            | Number of atoms in the BasisSet |
|`nbas`       | `Int32`                            | Number of basis functions (Note, not equal the number of BasisFunction objects) |
|`nshells`    | `Int32`                            | Number of shells, i.e. BasisFunction objects |
|`lc_atoms`   | `Array{Int32,1}`                   | Integer array maping data to libcint |
|`lc_bas`     | `Array{Int32,1}`                   | Integer array maping data to libcint |
|`lc_env`     | `Array{Float64,1}`                 | Float64 array maping data to libcint |

# Example

Build a basis set from default options
```julia
julia> bset = Fermi.GaussianBasis.BasisSet()
sto-3g Basis Set
Number of shells: 5
Number of basis:  7

O: 1s 2s 1p 
H: 1s 
H: 1s
```
The BasisSet object can be accessed as two-dimensional array.
```julia
julia> bset[2] # Show all basis for the second atom (H 1s)
1-element Vector{Fermi.GaussianBasis.BasisFunction}:
 S shell with 1 basis built from 3 primitive gaussians

χ₀₀  =    0.9817067283⋅Y₀₀⋅exp(-3.425250914⋅r²)
     +    0.9494640079⋅Y₀₀⋅exp(-0.6239137298⋅r²)
     +    0.2959064597⋅Y₀₀⋅exp(-0.168855404⋅r²)
julia> bset[1,2] # Show the second basis for the first atom (O 2s)
S shell with 1 basis built from 3 primitive gaussians

χ₀₀  =    0.8486970052⋅Y₀₀⋅exp(-5.033151319⋅r²)
     +    1.1352008076⋅Y₀₀⋅exp(-1.169596125⋅r²)
     +    0.8567529838⋅Y₀₀⋅exp(-0.38038896⋅r²)
```
You can also create your own crazy mix!
```julia
julia> @molecule {
    H 0.0 0.0 0.0
    H 0.0 0.0 0.7
} 
julia> mol = Fermi.Molecule()
# Lets create an S basis function for H
julia> s = Fermi.GaussianBasis.BasisFunction(0, [0.5215367271], [0.122])
# Now a P basis function
julia> p = Fermi.GaussianBasis.BasisFunction(1, [1.9584045349], [0.727])
# Now to create the basis set we need a dictionary maping atoms do basis functions
# Let's fetch the atom list from the molecule object
julia> Atoms = mol.atoms
# Now create the desired mapping as a dictionary. We want the first hydrogen
# to use only the S function, whereas the second will use both S and P
julia> shells = Dict(
    Atoms[1] => [s],
    Atoms[2] => [s,p]
)
# Now we can create the basis set!
julia> Fermi.GaussianBasis.BasisSet(mol, "My Crazy", shells)
My Crazy Basis Set
Number of shells: 3
Number of basis:  5

H: 1s 
H: 1s 1p
```
"""
struct BasisSet
    molecule::Molecule
    basis_name::String
    basis::Dict{Atom,Array{BasisFunction,1}}
    natoms::Cint
    nbas::Cint
    nshells::Cint
    lc_atoms::Array{Cint,1}
    lc_bas::Array{Cint, 1}
    lc_env::Array{Cdouble,1}
end

function BasisSet()
    mol = Molecule()
    bname = Options.get("basis")
    BasisSet(mol, bname)
end

function BasisSet(mol::Molecule, basis_name::String)

    output("\n  ⇒ Preparing new Basis Set: {:s}", basis_name)
    shells = Dict{Atom, Array{BasisFunction,1}}()
    for A in mol.atoms
        shells[A] = read_basisset(basis_name, Molecules.symbol(A))
    end
    BasisSet(mol, basis_name, shells)
end

function BasisSet(mol::Molecule, basis_name::String, shells::Dict{Atom, Array{BasisFunction,1}})

    atoms = mol.atoms
    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = length(atoms)
    nbas = 0
    nshells = 0
    nexps = 0
    nprims = 0

    for A in atoms
        basis = shells[A]
        for b in basis
            nshells += 1
            nbas += 2*b.l + 1
            nexps += length(b.exp)
            nprims += length(b.coef)
        end
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

"""
    Fermi.GaussianBasis.gto_norm(n::Signed, a::AbstractFloat)

Function that returns the normalization factor for a gaussian basis function rⁿ⋅exp(-a⋅r²)
"""
function gto_norm(n::Signed, a::AbstractFloat)
   # normalization factor of function rⁿ exp(-ar²)
    s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
    return √s
end

"""
    Fermi.GaussianBasis.normalize_basisfunction!(B::BasisFunction)

Modify the BasisFunction object by normalizing each primitive.
"""
function normalize_basisfunction!(B::BasisFunction)
    for i = eachindex(B.coef)
        B.coef[i] *= gto_norm(B.l, B.exp[i])
    end
end

function getindex(B::BasisSet, N::Int)
    A = B.molecule.atoms[N]
    return B.basis[A]
end

function getindex(B::BasisSet, i::Int, j::Int)
    A = B.molecule.atoms[i]
    return B.basis[A][j]
end

function string_repr(B::BasisFunction)
    # Generate Unicode symbol for sub number
    l_sub = Char(0x2080 + B.l)

    # Unicode for superscript is a bit messier, so gotta use control flow
    l_sup = B.l == 1 ? Char(0x00B9) :
            B.l in [2,3] ? Char(0x00B1 + B.l) :
            Char(0x2070 + B.l)

    nbas = 2*B.l + 1
    mvals = collect(-B.l:B.l)
    nprim = length(B.exp)

    # Reverse Dict(symbol=>num) to get Symbols from B.l
    Lsymbol = Dict(value => key for (key, value) in Fermi.GaussianBasis.AMDict)[B.l]
    out = "$(Lsymbol) shell with $nbas basis built from $nprim primitive gaussians\n\n"
    for m in mvals
        # Add sub minus sign (0x208B) if necessary
        m_sub = m < 0 ? Char(0x208B)*Char(0x2080 - m) : Char(0x2080 + m)
        out *= format("{:<4s} = ","χ$(l_sub)$m_sub")
        for i in eachindex(B.coef)

            if i > 1
                out *= B.coef[i] > 0 ? "\n     + " : "\n     - "
            end

            #out *= "$(abs(B.coef[i]))⋅Y$(l_sub)$m_sub"
            out *= format("{:>15.10f}⋅Y$(l_sub)$m_sub", abs(B.coef[i]))

            if B.l != 0 
                out *= "⋅r$l_sup"
            end

            out *= "⋅exp(-$(B.exp[i])⋅r²)"
        end
        out *="\n\n"
    end
    return strip(out)
end

function string_repr(B::BasisSet)
    out  =  "$(B.basis_name) Basis Set\n"
    out *= "Number of shells: $(B.nshells)\n"
    out *= "Number of basis:  $(B.nbas)\n\n"

    l_to_symbol = Dict(
        0 => "s",
        1 => "p",
        2 => "d",
        3 => "f",
        4 => "g",
        5 => "h",
        6 => "i",
    )

    for A in B.molecule.atoms
        # Count how many times s,p,d appears for numbering
        count = zeros(Int16, 7)
        out *= "$(Molecules.symbol(A)): "
        for bfs in B.basis[A]
            L = bfs.l
            count[L+1] += 1
            out *= "$(count[L+1])$(l_to_symbol[L]) "
        end
        out *="\n"
    end

    return strip(out)
end

# Pretty printing
function show(io::IO, ::MIME"text/plain", X::T) where T<:Union{BasisFunction, BasisSet}
    print(io, string_repr(X))
end

end #module