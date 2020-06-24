abstract type AbstractMolecule end
begin
    mutable struct Molecule <: AbstractMolecule
        lmolecule::M where M <: Lints.Molecule
        atoms::Array{Fermi.Atom,1} 
        charge::Int8
        multiplicity::Int8
    end
    struct NullMolecule <: AbstractMolecule end
end

function Molecule(atoms::Array{Fermi.Atom},charge::Int8)
    # guess low spin multiplicity
    Molecule(atoms,charge,mult)
end
function Molecule(atoms::Array{Fermi.Atom},charge::Int8,mult::Int8)
    # make xyz string
    # write xyz string to /tmp/fermi_<UUID>.xyz
    # load Lints.Molecule from ^^
    Molecule(lmol,atoms,charge,mult)
end
