"""
    Fermi.Integrals

Module to compute integrals using Lints.jl
"""
module Integrals

using Fermi
using Fermi.Error
using Fermi.Options
using Fermi.Geometry: Molecule
using LinearAlgebra
using TensorOperations

import Base: getindex, setindex!, delete!

export IntegralHelper
export delete!

include("Lints.jl")

"""
    IntegralHelper{T}

Structure to assist with computing and storing integrals. 
Accesss like a dictionary e.g.,
    ints["S"]

A key is associated with each type of integral

    "S"           -> AO overlap integral
    "T"           -> AO electron kinetic energy integral
    "V"           -> AO electron-nuclei attraction integral
    "ERI"         -> AO electron repulsion integral
    "JKERI"       -> AO JK density fitted electron repulsion integral
    "RIERI"       -> AO RI density fitted electron repulsion integral

# Fields
    mol                         Associated Fermi Molecule object
    basis                       Basis set used within the helper
    aux                         Auxiliar basis set used in density fitting
    cache                       Holds integrals already computed 
    normalize                   Do normalize integrals? `true` or `false`
"""
struct IntegralHelper{T<:AbstractFloat,E<:AbstractERI,O<:AbstractOrbitals}
    molecule::Molecule
    orbitals::O
    basis::String
    cache::Dict{String,FermiMDArray{T}} 
    eri_type::E
    normalize::Bool
end

function IntegralHelper(;molecule = Molecule(), orbitals = AtomicOrbitals(), 
                           basis = Options.get("basis"), normalize = false)

    # Check if density-fitting is requested
    if Options.get("df")
        # If the associated orbitals are AtomicOrbitals and DF is requested, JKFIT is set by default
        # Otherwise, the ERI type will be RIFIT
        eri_type = orbitals === AtomicOrbitals() ? JKFIT() : RIFIT()
    else
        eri_type = Chonky()
    end

    # Starts an empty cache

    precision = Options.get("precision")
    if precision == "single"
        cache = Dict{String, FermiMDArray{Float32}}() 
    elseif precision == "double"
        cache = Dict{String, FermiMDArray{Float64}}() 
    else
        throw(InvalidFermiOption("precision can only be `single` or `double`. Got $precision"))
    end

    # Return IntegralHelper object
    IntegralHelper(molecule, orbitals, basis, cache, eri_type, normalize)
end

function IntegralHelper(I::IntegralHelper{T,E,O}, A::AtomicOrbitals) where {T<:AbstractFloat,E<:AbstractERI,O<:AbstractOrbitals}
    cache = Dict{String, FermiMDArray{T}}() 
    return IntegralHelper{T,E,AtomicOrbitals}(I.molecule, A, I.basis, cache, E(), I.normalize)
end

# Clears cache and change normalize key
function normalize!(I::IntegralHelper,normalize::Bool)
    if I.normalize != normalize
        I.normalize = normalize
        for entry in keys(I.cache)
            delete!(I.cache, entry)
        end
    end
end

"""
    getindex(I::IntegralHelper,entry::String)

Called when `helper["foo"]` syntax is used. If the requested entry already
exists, simply return the entry. If not, compute the requested entry.
"""
function getindex(I::IntegralHelper,entry::String)
    if haskey(I.cache, entry)
        return I.cache[entry]
    else
        compute!(I, entry)
        return I.cache[entry]
    end
end

function setindex!(I::IntegralHelper, A::FermiMDArray, key::String)
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

include("AtomicIntegrals.jl")
include("ROIntegrals.jl")

end # Module