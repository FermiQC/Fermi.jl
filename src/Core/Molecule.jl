using Molecules
import Molecules: Atom

export Molecule, Atom

using Fermi.Options
using Fermi.PhysicalConstants: atomic_number, bohr_to_angstrom

using LinearAlgebra
using Formatting

import Base: show

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
julia> Molecule()
Molecule:

O    1.209153654800    1.766411818900   -0.017161397200
H    2.198480007500    1.797710062700    0.012116171900
H    0.919788188200    2.458018557000    0.629793883200


Charge: 0   Multiplicity: 1   
Nuclear repulsion:    8.8880641737

julia> Molecule(charge=2, multiplicity=3)
Molecule:

O    1.209153654800    1.766411818900   -0.017161397200
H    2.198480007500    1.797710062700    0.012116171900
H    0.919788188200    2.458018557000    0.629793883200


Charge: 2   Multiplicity: 3   
Nuclear repulsion:    8.8880641737
```
"""
struct Molecule
    atoms::Vector{Atom}
    charge::Int
    multiplicity::Int
    Nα::Int
    Nβ::Int
end

function Molecule(;
    molstring = Options.get("molstring"),
    unit = Symbol(Options.get("unit")),
    charge = Options.get("charge"),
    multiplicity = Options.get("multiplicity")
    )

    atoms = Molecules.parse_string(molstring, unit=unit)
    Molecule(atoms, charge, multiplicity)
end

function Molecule(atoms::Vector{T}, charge::Int, multiplicity::Int) where T <: Atom

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
    out =  Molecule(atoms, charge, multiplicity, Nα, Nβ)
    return out
end

"""
    Fermi.string_repr(M::Molecule)

Returns a nicely formatted string with all the molecule's information
"""
function string_repr(M::Molecule)
    out = ""
    out = out*format("Molecule:\n\n")
    out = out*format(Molecules.get_xyz(M.atoms))
    out = out*format("\n")
    out = out*format("\nCharge: {}   ", M.charge)
    out = out*format("Multiplicity: {}   \n", M.multiplicity)
    return out
end

# Pretty printing
function show(io::IO, ::MIME"text/plain", X::Molecule)
    print(io, string_repr(X))
end