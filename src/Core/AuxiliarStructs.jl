using Fermi.Geometry: Molecule

export AbstractOrbitals, AtomicOrbitals, AbstractRestrictedOrbitals, AbstractUnrestrictedOrbitals, GeneralRestrictedOrbitals
export AbstractWavefunction
export AbstractERI, AbstractDFERI, JKFIT, RIFIT, Chonky

"""
    Fermi.AbstractWavefunction

Abstract type common to all wave functions.

_struct tree:_

**AbstractWavefunction**  (Top level)
"""
abstract type AbstractWavefunction end

abstract type AbstractERI end
abstract type AbstractDFERI <: AbstractERI end

struct JKFIT <: AbstractDFERI end
struct RIFIT <: AbstractDFERI end
struct Chonky <:AbstractERI end
