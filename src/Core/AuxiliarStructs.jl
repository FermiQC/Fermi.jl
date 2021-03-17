using Fermi.Geometry: Molecule

"""
    Fermi.AbstractWavefunction

Abstract type common to all wave functions.

_struct tree:_

**AbstractWavefunction**  (Top level)
"""
abstract type AbstractWavefunction end

"""
    Fermi.AbstractOrbitals

Abstract type common to all orbitals

_struct tree:_

**AbstractOrbitals**  (Top level)
"""
abstract type AbstractOrbitals end

"""
    Fermi.GeneralRestrictedOrbitals

Dummy orbital struct that can be used to hold any type of Restricted orbitals.
Can be used to read in custom orbitals as

    Orbs = GeneralRestrictedOrbitals(X)

where X is an Array object. The molecule and basis will be deduced from the current options.
Alternatively, one can pass these informations explicitly:

    Orbs = GeneralRestrictedOrbitals(X, molecule=mol, name="myorbitals", basis="cc-pvdz")

`name` and `basis` are Strings, whereas `mol` is a `Fermi.Geometry.Molecule` object.

# Fields:

    name      String with a label for the orbital
    basis     String indicating the basis set used to construct the orbitals
    aux       String indicating the auxiliar basis set used for density fitting
    molecule  Molecule object for which the orbitals were constructed
    C         NxN AbstractArray with the AO(lines) → MO(orbitals) coefficients

_struct tree:_

**GeneralRestrictedOrbitals** <: AbstractOrbitals
"""
struct GeneralRestrictedOrbitals{T} <: AbstractOrbitals 
    name::String
    basis::String
    aux::String
    molecule::Molecule
    C::AbstractArray{T,2}
end

function GeneralRestrictedOrbitals(C::AbstractArray{T,2}; mol=nothing, name="Custom", basis="undef", aux="undef") where T <: AbstractFloat

    mol === nothing ? mol = Fermi.Geometry.Molecule() : nothing
    basis == "undef" ? basis = Fermi.Options.get("basis") : nothing
    aux == "undef" ? basis = Fermi.Options.get("aux") : nothing

    GeneralRestrictedOrbitals{T}(name, basis, mol, C)
end
