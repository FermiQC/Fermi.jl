"""
    Fermi.AbstractWavefunction

Abstract type common to all wave functions.

_struct tree:_

**AbstractWavefunction**  (Top level)
"""
abstract type AbstractWavefunction end

"""
    Fermi.AbstractReferenceWavefunction

Abstract type common to all reference wave functions.

_struct tree:_

**AbstractReferenceWavefunction** <: AbstractWavefunction
"""
abstract type AbstractReferenceWavefunction <: AbstractWavefunction end

"""
    Fermi.AbstractCorrelatedWavefunction

Abstract type common to all correlated wave functions.

_struct tree:_

**AbstractCorrelatedWavefunction** <: AbstractWavefunction
"""
abstract type AbstractCorrelatedWavefunction <: AbstractWavefunction end

