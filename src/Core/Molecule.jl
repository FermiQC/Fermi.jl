using Molecules
import Molecules: Atom, Molecule

export Molecule, Atom

using Fermi.Options
using Fermi.PhysicalConstants: atomic_number, bohr_to_angstrom

using LinearAlgebra
using Formatting

import Base: show

"""
    Fermi.Molecule

Construct an object storing information about a molecule (group of atoms).
Note that this object is imported from the `Molecules` package, which is a dependency of Fermi.

# Fields
    
    atoms         Array with Fermi.Atom objects
    charge        Charge of the molecule
    multiplicity  Multiplicity ``(2Mₛ + 1)``
    Nα            Number of Alpha electrons
    Nβ            Number of Beta electrons

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

julia> Molecule(charge=2, multiplicity=3)
Molecule:

O    1.209153654800    1.766411818900   -0.017161397200
H    2.198480007500    1.797710062700    0.012116171900
H    0.919788188200    2.458018557000    0.629793883200


Charge: 2   Multiplicity: 3   
```
"""
function Molecule(;
    molstring = Options.get("molstring"),
    unit = Symbol(Options.get("unit")),
    charge = Options.get("charge"),
    multiplicity = Options.get("multiplicity")
    )

    atoms = Molecules.parse_string(molstring, unit=unit)
    Molecule(atoms, charge, multiplicity)
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