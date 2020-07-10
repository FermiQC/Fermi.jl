"""
    Fermi.MolHelper

Module handling molecule and atoms data

# exports:

    Atom      Object storing information about an atom
    Molecule  Object storing information about a molecule (group of atoms)
"""
module MolHelper

export Molecule
export Atom

using Fermi
using Fermi.PhysicalConstants: atomic_number, bohr_to_angstrom, angstrom_to_bohr
using LinearAlgebra

"""
    Fermi.Atom

Object storing information about an atom.

# Fields:
    
    AtomicSymbol    Atomic symbol
    Z               Atomic number
    xyz             xyz array in Angstrom
"""
struct Atom
    AtomicSymbol::String
    Z::Int
    xyz::Array{Float64,1}
end

"""
    Fermi.Molecule

Object storing information about a molecule (group of atoms).

# Fields:
    
    atoms        Array with Fermi.Atom objects
    charge       Charge of the molecule
    multiplicity Multiplicity ``(2Ms + 1)``
    Nα           Number of Alpha electrons
    Nβ           Number of Beta electrons
    Vnuc         Nuclear repulsion energy
"""
struct Molecule
    atoms::Array{Atom,1}
    charge::Int
    multiplicity::Int
    Nα::Int
    Nβ::Int
    Vnuc::Float64
end

"""
    Fermi.Molecule()

Uses data stored in Fermi.CurrentOptions to return a Molecule object.

"""
function Molecule()
    molstring = Fermi.CurrentOptions["molstring"]
    unit = Fermi.CurrentOptions["unit"]
    Molecule(molstring, unit)
end

"""
    Fermi.Molecule(molstring::String, unit::String="angstrom")

Uses data stored in Fermi.CurrentOptions for charge and multiplicity to build
a Molecule object given XYZ string.
"""
function Molecule(molstring::String, unit::String="angstrom")
    charge = Fermi.CurrentOptions["charge"]
    multiplicity = Fermi.CurrentOptions["multiplicity"]
    Molecule(molstring, charge, multiplicity, unit)
end

"""
    Fermi.Molecule(molstring::String, charge::Int, multiplicity::Int, unit::String="angstrom")

Produces a Molecule object from the XYZ string.
"""
function Molecule(molstring::String, charge::Int, multiplicity::Int, unit::String="angstrom")
    atoms = get_atoms(molstring, unit=unit)
    Molecule(atoms, charge, multiplicity)
end

"""
    Fermi.Molecule(molstring::String, charge::Int, multiplicity::Int, unit::String="angstrom")

Produces a Molecule object from an array of Atom objects.
"""
function Molecule(atoms::Array{Atom,1}, charge::Int, multiplicity::Int)
    
    # Compute Nuclear repulsion
    Vnuc = 0.0
    for i in eachindex(atoms)
        for j in 1:(i-1)
            @inbounds Vnuc += nuclear_repulsion(atoms[i], atoms[j])
        end
    end

    # Compute number of electrons

    nelec = 0
    for i in eachindex(atoms)
        @inbounds nelec += atoms[i].Z
    end

    nelec -= charge

    if nelec ≤ 0
        throw(Fermi.InvalidFermiOption("Invalid charge ($charge) for given molecule"))
    end

    αexcess = multiplicity-1

    if isodd(nelec) != isodd(αexcess)
        throw(Fermi.InvalidFermiOption("Incompatible charge $(charge) and multiplicity $(multiplicity)"))
    end

    Nβ = (nelec - αexcess)/2
    Nα = nelec - Nβ
    
    return Molecule(atoms, charge, multiplicity, Nα, Nβ, Vnuc)
end

"""
    Fermi.MolHelper.nuclear_repulsion(A::Atom, B::Atom)

Returns the repulsion energy between the two given atoms.
"""
function nuclear_repulsion(A::Atom, B::Atom)
    return (A.Z*B.Z)/(√((A.xyz.-B.xyz)⋅(A.xyz.-B.xyz))*angstrom_to_bohr)
end

"""
    Fermi.get_atoms(molstring::String; unit::String="angstrom")

From a XYZ string, produces an array of Atom objects.
"""
function get_atoms(molstring::String; unit::String="angstrom")
    
    if unit == "bohr"
        conv = bohr_to_angstrom
    elseif unit == "angstrom"
        conv = 1.0
    else
        error("Invalid unit given to Fermi.Molecule.get_atoms: $unit")
    end

    atoms = Atom[]
    exp = r"(\w{1,2})\s+([-+]??\d+\.\d+\s+[-+]??\d+\.\d+\s+[-+]??\d+\.\d+)"
    for line in split(molstring, "\n")
        m = match(exp, line)
        if m != nothing
            # Convert SubString to String
            AtomicSymbol = String(m.captures[1])
            Z = atomic_number(AtomicSymbol)

            # Convert Substring to String and then String to Float64 Array
            xyz = parse.(Float64, split(String(m.captures[2]))).*conv

            push!(atoms, Atom(AtomicSymbol, Z, xyz))
        end
    end
    return atoms
end

"""
    Fermi.get_xyz(M::Molecule)

Returns a XYZ string in angstrom for the given Molecule.
"""
function get_xyz(M::Molecule)

    molstring = ""

    for i in eachindex(M.atoms)

        A = M.atoms[i]
        xyz = ["$(A.xyz[1])", "$(A.xyz[2])", "$(A.xyz[3])"]

        # Format the string because of OCD

        ## Add space in front
        for i in eachindex(xyz)
            if !(xyz[i][1] in "+-")
                xyz[i] = " "*xyz[i]
            end
        end

        ## Equal the length
        L = max(map(length, xyz)...)

        for i in eachindex(xyz)
            while length(xyz[i]) < L
                xyz[i] = xyz[i]*"0"
            end
        end
        molstring *= A.AtomicSymbol*"   "*join(xyz, "   ")*"\n"
    end

    return String(strip(molstring))
end

end #Module
