export JKFIT, RIFIT, Chonky, AbstractDFERI, AbstractERI, SparseERI

"""
    Fermi.Integrals.AbstractERI

Abstract type common to all ERI types
"""
abstract type AbstractERI end

"""
    Fermi.Integrals.AbstractDFERI

Abstract type common to all density-fitted ERI types

# Struct tree

**AbstractDFERI** <: AbstractERI
"""
abstract type AbstractDFERI <: AbstractERI end

"""
    Fermi.Integrals.JKFIT

Concrete type representing a density-fitted ERI using a JKFIT auxiliar basis set.

# Fields

| Name       | Type       | Description |
|:-----------|:-----------|:----------------------------------------------------|
|`basisset`  | `BasisSet` | BasisSet object for the auxiliar JK basis |

# Examples

The JKFIT structure is used to build a IntegralHelper
```julia
julia> jk = Fermi.Integrals.JKFIT()
julia> ints = Fermi.Integrals.IntegralHelper(eri_type=jk) 
```

# Struct tree

JKFIT <: AbstractDFERI <: AbstractERI
"""
struct JKFIT <: AbstractDFERI
    basisset::BasisSet
end

function JKFIT(mol::Molecule = Molecule())

    auxjk = Options.get("jkfit")
    # If aux is auto, determine the aux basis from the basis
    if auxjk == "auto"
        basis = Options.get("basis")
        dunning_name = Regex("cc-pv.z")
        auxjk = occursin(dunning_name, basis) ? basis*"-jkfit" : "cc-pvqz-jkfit"
    end

    return JKFIT(mol, auxjk)
end

function JKFIT(mol::Molecule, basis::String)
    return JKFIT(BasisSet(basis, mol.atoms))
end

"""
    Fermi.Integrals.RIFIT

Concrete type representing a density-fitted ERI using a RIFIT auxiliar basis set.

# Fields

| Name       | Type       | Description |
|:-----------|:-----------|:----------------------------------------------------|
|`basisset`  | `BasisSet` | BasisSet object for the auxiliar RI basis |

# Struct tree

RIFIT <: AbstractDFERI <: AbstractERI
"""
struct RIFIT <: AbstractDFERI 
    basisset::BasisSet
end

function RIFIT(mol::Molecule = Molecule())

    auxri = Options.get("rifit")
    # If aux is auto, determine the aux basis from the basis
    if auxri == "auto"
        basis = Options.get("basis")
        dunning_name = Regex("cc-pv.z")
        auxri = occursin(dunning_name, basis) ? basis*"-rifit" : "cc-pvqz-rifit"
    end

    return RIFIT(mol, auxri)
end

function RIFIT(mol::Molecule, basis::String)
    return RIFIT(BasisSet(basis, mol.atoms))
end

struct Chonky <:AbstractERI end
struct SparseERI <: AbstractERI end