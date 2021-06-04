module Orbitals

using Fermi
using Fermi.Geometry
using Fermi.GaussianBasis
using TensorOperations

export AbstractOrbitals, AtomicOrbitals, AbstractRestrictedOrbitals, AbstractUnrestrictedOrbitals
export GeneralRestrictedOrbitals
export RHFOrbitals

"""
    Fermi.AbstractOrbitals

Abstract type common to all orbitals

_struct tree:_

**AbstractOrbitals**  (Top level)
"""
abstract type AbstractOrbitals end

struct AtomicOrbitals <: AbstractOrbitals 
    basisset::BasisSet
end

function AtomicOrbitals()
    return AtomicOrbitals(BasisSet())
end

function AtomicOrbitals(mol:: Molecule, basis::String)
    return AtomicOrbitals(BasisSet(mol, basis))
end

abstract type AbstractRestrictedOrbitals <: AbstractOrbitals end
abstract type AbstractUnrestrictedOrbitals <: AbstractOrbitals end

"""
    Fermi.GeneralRestrictedOrbitals

Dummy orbital struct that can be used to hold any type of Restricted orbitals.
Can be used to read in custom orbitals as

    Orbs = GeneralRestrictedOrbitals(X)

where X is an Array object. The molecule and basis will be deduced from the current options.
Alternatively, one can pass these informations explicitly:

    Orbs = GeneralRestrictedOrbitals(X, molecule=mol, name="myorbitals", basis="cc-pvdz")

`name` and `basis` are Strings, whereas `mol` is a `Fermi.Geometry.Molecule` object.

# Fields

    name      String with a label for the orbital
    basis     String indicating the basis set used to construct the orbitals
    molecule  Molecule object for which the orbitals were constructed
    C         NxN AbstractArray with the AO(lines) → MO(orbitals) coefficients

_struct tree:_

**GeneralRestrictedOrbitals** <: AbstractOrbitals
"""
struct GeneralRestrictedOrbitals{T} <: AbstractRestrictedOrbitals 
    molecule::Molecule
    basis::String
    sd_energy::T
    C::AbstractArray{T,2}
end

function GeneralRestrictedOrbitals(C::AbstractArray{T,2}; mol=nothing, basis="undef", sd_energy=zero(T)) where T <: AbstractFloat

    mol === nothing ? mol = Fermi.Geometry.Molecule() : nothing
    basis == "undef" ? basis = Fermi.Options.get("basis") : nothing

    GeneralRestrictedOrbitals{T}(mol, basis, sd_energy, C)
end

"""
    Fermi.HartreeFock.RHFOrbitals

Struct holding information about Restricted Hartree--Fock orbitals

# Fields

    molecule   Molecule object associated with the orbitals
    basis      Basis set used to compute the orbitals
    eps        Orbital energies, i.e. diagonal of the Fock matrix
    C          Coefficients of the AO->MO transformation matrix
"""
struct RHFOrbitals <: AbstractRestrictedOrbitals
    molecule::Molecule
    basis::String
    eps::AbstractArray{Float64,1}
    sd_energy::Float64
    C::AbstractArray{Float64,2}
end

end # Module