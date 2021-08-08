"""
    Fermi.Geometry

Module handling molecule and atoms data

# Objects

    Atom      Object storing information about an atom
    Molecule  Object storing information about a molecule (group of atoms)
"""
module Geometry

export Molecule
export Atom

using Fermi
using Fermi.Options
using Fermi.PhysicalConstants: atomic_number, bohr_to_angstrom

using LinearAlgebra
using Formatting

import Base: show

"""
    Fermi.Atom

Object storing information about an atom.

# Fields

| Name  | Description |
|:------|:-----------------------------------------------------------|
|`AtomicSymbol` |  Atomic symbol          |
|`Z`            |  Atomic number          |
|`xyz`          |  xyz array in Angstrom  |
"""
struct Atom
    AtomicSymbol::String
    Z::Int
    xyz::Tuple{Float64,Float64,Float64}
end

"""
    Fermi.Molecule

Object storing information about a molecule (group of atoms).

# Fields
    
    atoms         Array with Fermi.Atom objects
    charge        Charge of the molecule
    multiplicity  Multiplicity ``(2Mₛ + 1)``
    Nα            Number of Alpha electrons
    Nβ            Number of Beta electrons
    Vnuc          Nuclear repulsion energy


# Examples:

A Molecule object can be created by providing the desired keyworded arguments: 

    molstring:      A string representing the XYZ of the molecule
    unit:           Unit used for the distance between atoms (Bohr or Angstrom)
    charge:         Molecular charge
    multiplicity    Molecular multiplicity

Any argument not given explicitly will be read from the Options.
```
julia> Fermi.Geometry.Molecule()
Molecule:

O    1.209153654800    1.766411818900   -0.017161397200
H    2.198480007500    1.797710062700    0.012116171900
H    0.919788188200    2.458018557000    0.629793883200


Charge: 0   Multiplicity: 1   
Nuclear repulsion:    8.8880641737

julia> Fermi.Geometry.Molecule(charge=2, multiplicity=3)
Molecule:

O    1.209153654800    1.766411818900   -0.017161397200
H    2.198480007500    1.797710062700    0.012116171900
H    0.919788188200    2.458018557000    0.629793883200


Charge: 2   Multiplicity: 3   
Nuclear repulsion:    8.8880641737
```
"""
struct Molecule
    atoms::Tuple
    charge::Int
    multiplicity::Int
    Nα::Int
    Nβ::Int
    Vnuc::Float64
end

function Molecule(;
    molstring = Options.get("molstring"),
    unit = Options.get("unit"),
    charge = Options.get("charge"),
    multiplicity = Options.get("multiplicity")
    )

    atoms = get_atoms(molstring, unit=unit)
    Molecule(atoms, charge, multiplicity)
end

function Molecule(atoms::Array{Atom,1}, charge::Int, multiplicity::Int)
    
    # Compute Nuclear repulsion
    Vnuc = 0.0
    for i in eachindex(atoms)
        for j in 1:(i-1)
            Vnuc += nuclear_repulsion(atoms[i], atoms[j])
        end
    end

    # Compute number of electrons
    nelec = -charge
    for i in eachindex(atoms)
        nelec += atoms[i].Z
    end

    # If the number of electrons turns out negative returns an error
    if nelec ≤ 0
        throw(FermiException("Invalid charge ($charge) for given molecule"))
    end

    # Mult = 2Ms + 1 thus the number of unpaired electrons (taken as α) is Mult-1 = 2Ms
    αexcess = multiplicity-1

    # The number of β electrons must be equal the number of doubly occupied orbitals (Ndo)
    # Ndo = (nelec - n_of_unpaired_elec)/2 this must be an integer
    if isodd(nelec) != isodd(αexcess)
        throw(FermiException("Incompatible charge $(charge) and multiplicity $(multiplicity)"))
    end

    Nβ = (nelec - αexcess)/2
    Nα = nelec - Nβ
    out =  Molecule(Tuple(atoms), charge, multiplicity, Nα, Nβ, Vnuc)
    return out
end

"""
    Fermi.get_atoms(molstring::String; unit::String="angstrom")

From a XYZ string, produces an array of Atom objects.
"""
function get_atoms(molstring::String; unit::String="angstrom")
    
    # Get a list of Atom objects from String
    if unit == "bohr"
        conv = bohr_to_angstrom
    elseif unit == "angstrom"
        conv = 1.0
    else
        throw(FermiException("unknown unit in molecule construction: $unit"))
    end

    atoms = Atom[]
    for line in split(strip(molstring), "\n")
        m = split(line)
        if length(m) != 4
            throw(FermiException("4 columns expected on XYZ string. Got $(length(m))"))
        end

        if m !== nothing
            # Convert SubString to String
            AtomicSymbol = String(m[1])
            Z = atomic_number(AtomicSymbol)

            # Convert Substring to String and then String to Float64 Array
            xyz = [0.0, 0.0, 0.0]
            try
                xyz .= parse.(Float64, m[2:4]).*conv
            catch ArgumentError
                throw(FermiException("Failed to process XYZ string: $(m[2:4])"))
            end

            push!(atoms, Atom(AtomicSymbol, Z, Tuple(xyz)))
        end
    end
    return atoms
end

"""
    Fermi.Geometry.nuclear_repulsion(A::Atom, B::Atom)

Returns the repulsion energy between atoms A and B.
"""
function nuclear_repulsion(A::Atom, B::Atom)
    return (A.Z*B.Z)/(√((A.xyz.-B.xyz)⋅(A.xyz.-B.xyz))/bohr_to_angstrom)
end

"""
    Fermi.get_xyz(M::Molecule)

Returns a XYZ string in angstrom for the given Molecule.
"""
function get_xyz(M::Molecule)

    molstring = ""

    for i in eachindex(M.atoms)

        A = M.atoms[i]
        molstring *= format("{}   {: 15.12f}   {: 15.12f}   {: 15.12f}\n", A.AtomicSymbol, A.xyz...)
    end
    return molstring
end

"""
    Fermi.string_repr(M::Molecule)

Returns a nicely formatted string with all the molecule's information
"""
function string_repr(M::Molecule)
    out = ""
    out = out*format("Molecule:\n\n")
    out = out*format(Fermi.Geometry.get_xyz(M))
    out = out*format("\n")
    out = out*format("\nCharge: {}   ", M.charge)
    out = out*format("Multiplicity: {}   \n", M.multiplicity)
    out = out*format("Nuclear repulsion: {:15.10f}", M.Vnuc)
    return out
end

"""
    Fermi.string_repr(A::Atom)

Returns a nicely formatted string for Atom object.
"""
function string_repr(A::Atom)
    out = "Symbol:           $(A.AtomicSymbol)\n"
    out = out*"Atomic Number:    $(A.Z)\n"
    out = out*"Position: $(format("{:2.5f}  {:2.5f}  {:2.5f}", A.xyz...))"
end

# Pretty printing
function show(io::IO, ::MIME"text/plain", X::T) where T<:Union{Molecule, Atom}
    print(io, string_repr(X))
end
end #Module