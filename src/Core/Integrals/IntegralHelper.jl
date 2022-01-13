"""
    Fermi.Integrals

Module to compute and manage molecular integrals.
"""
module Integrals

using Fermi
using Fermi.Options
using Fermi.Orbitals
using LinearAlgebra
using GaussianBasis
using TensorOperations

import Fermi: string_repr
import Base: getindex, setindex!, delete!, show
import GaussianBasis: BasisSet

export IntegralHelper, delete!, mo_from_ao!
export BasisSet, BasisFunction

# Expand BasisSet methods to handle molecules.
BasisSet(name::String, mol::Molecule) = BasisSet(name, mol.atoms)

include("ERITypes.jl")

"""
    IntegralHelper

Manager for integrals computation and storage.
Accesss like a dictionary e.g.,
```
julia> ints = Fermi.Integrals.IntegralHelper()
julia> ints["S"] #Returns overlap matrix
```
A key is associated with each type of integral

    "S"           -> Overlap integral
    "T"           -> Electron kinetic energy integral
    "V"           -> Electron-nuclei attraction integral
    "ERI"         -> Electron repulsion integral

# Fields
    molecule                    Molecule object
    orbitals                    Orbitals used in the integral computation
    basis                       Basis set name
    cache                       Holds computed integrals 
    eri_type                    Defines how electron repulsion integrals are handled 
"""
struct IntegralHelper{T<:AbstractFloat,E<:AbstractERI,O<:AbstractOrbitals}
    molecule::Molecule
    orbitals::O
    basis::String
    cache::Dict{String,AbstractArray{T}} 
    eri_type::E
end

function IntegralHelper(x...;k...)

    precision = Options.get("precision")
    if precision == "single"
        IntegralHelper{Float32}(x...; k...)
    elseif precision == "double"
        IntegralHelper{Float64}(x...; k...)
    else
        throw(FermiException("precision can only be `single` or `double`. Got $precision"))
    end
end

function IntegralHelper{T}(;molecule = Molecule(), orbitals = nothing, 
                           basis = Options.get("basis"), eri_type=nothing) where T<:AbstractFloat

    if orbitals === nothing
        bset = BasisSet(basis, molecule)
        orbitals = AtomicOrbitals(bset)
    end

    # If the associated orbitals are AtomicOrbitals and DF is requested, JKFIT is set by default
    if Options.get("df") && orbitals isa AtomicOrbitals && eri_type === nothing
        eri_type = JKFIT(molecule)

    # If df is requested, but the orbitals are not AtomicOrbitals then RIFIT is set by default
    elseif Options.get("df") && eri_type === nothing

        eri_type = RIFIT(molecule)

    # Else, if eri_type is not an AbstractERI object, Chonky is used. This will overwrite invalid entries for eri_type
    elseif !(typeof(eri_type) <: AbstractERI)
        eri_type = orbitals isa AtomicOrbitals ? SparseERI() : Chonky()
    end

    cache = Dict{String, AbstractArray{T}}() 
    IntegralHelper(molecule, orbitals, basis, cache, eri_type)
end

function IntegralHelper{T}(bset::BasisSet, eri_type=nothing) where T<:AbstractFloat
    IntegralHelper{T}(
        molecule = Molecule(bset.atoms, Options.get("charge"), Options.get("multiplicity")),
        orbitals = AtomicOrbitals(bset),
        basis = bset.name,
        eri_type=eri_type
    )
end

"""
    getindex(I::IntegralHelper,entry::String)

Called when `helper["foo"]` syntax is used. If the requested entry already
exists, simply return the entry. If not, compute the requested entry.
"""
function getindex(I::IntegralHelper, entry::String)
    if haskey(I.cache, entry)
        return I.cache[entry]
    else
        compute!(I, entry)
        return I.cache[entry]
    end
end

function setindex!(I::IntegralHelper, A::AbstractArray, key::String)
    I.cache[key] = A
end

function delete!(I::IntegralHelper, keys...)
    for k in keys
        delete!(I.cache, k)
    end
    GC.gc()
end

function delete!(I::IntegralHelper)
    delete!(I, keys(I.cache)...)
end

# Integrals specific to orbital types
include("AtomicIntegrals.jl")

# MO conversions for Restricted Orbitals
include("ROIntegrals/ROIntegrals.jl")

# Pretty printing
function string_repr(X::IntegralHelper{T,E,O}) where {T,E,O}
    eri_string = replace("$E", "Fermi.Integrals."=>"")
    orb_string = replace("$O", "Fermi.Orbitals."=>"")
    out = ""
    out = out*" ⇒ Fermi IntegralHelper\n"
    out = out*" ⋅ Data Type:                 $(T)\n"
    out = out*" ⋅ Basis:                     $(X.basis)\n"
    out = out*" ⋅ ERI:                       $(eri_string)\n"
    out = out*" ⋅ Orbitals:                  $(orb_string)\n"
    cache_str = ""
    for k in keys(X.cache)
        cache_str *= k*" "
    end
    out = out*" ⋅ Stored Integrals:          $(cache_str)"
    return out
end

function string_repr(X::AbstractDFERI)
    eri_string = replace("$(typeof(X))", "Fermi.Integrals."=>"")
    return "$(eri_string): $(X.basisset.name)"
end

function show(io::IO, ::MIME"text/plain", X::T) where T<:Union{IntegralHelper,JKFIT,RIFIT}
    print(io, string_repr(X))
end

end # Module