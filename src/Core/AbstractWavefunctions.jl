"""
    Fermi.AbstractWavefunction

Fermi abstract type common to all wave functions.

_struct tree:_

**AbstractWavefunction**  (Top level)
"""
abstract type AbstractWavefunction end

"""
    Fermi.AbstractReferenceWavefunction

Fermi abstract type common to all reference wave functions.

_struct tree:_

**AbstractReferenceWavefunction** <: AbstractWavefunction
"""
abstract type AbstractReferenceWavefunction <: AbstractWavefunction end

"""
    Fermi.AbstractCorrelatedWavefunction

Fermi abstract type common to all correlated wave functions.

_struct tree:_

**AbstractCorrelatedWavefunction** <: AbstractWavefunction
"""
abstract type AbstractCorrelatedWavefunction <: AbstractWavefunction end

