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
using Fermi.Orbitals

import Base: getindex, setindex!, delete!, show

export IntegralHelper, delete!, mo_from_ao!

include("../Backend/Lints.jl")

"""
    IntegralHelper{T}

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
    cache                       Holds integrals computed 
    eri_type                    Defines whether/how density-fitting is done
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

# Pretty printing
function string_repr(X::IntegralHelper{T,E,O}) where {T,E,O}
    out = ""
    out = out*" ⇒ Fermi IntegralHelper\n"
    out = out*" ⋅ Data Type:                 $(T)\n"
    out = out*" ⋅ Basis:                     $(X.basis)\n"
    out = out*" ⋅ ERI:                       $(E)\n"
    out = out*" ⋅ Orbitals:                  $(O)\n"
    cache_str = ""
    for k in keys(X.cache)
        cache_str *= k*" "
    end
    out = out*" ⋅ Stored Integrals:          $(cache_str)"
    return out
end

function show(io::IO, ::MIME"text/plain", X::IntegralHelper)
    print(string_repr(X))
end


function IntegralHelper(;molecule = Molecule(), orbitals = AtomicOrbitals(), 
                           basis = Options.get("basis"), normalize = false, eri_type = nothing)

    precision = Options.get("precision")
    if precision == "single"
        IntegralHelper{Float32}(molecule=molecule, orbitals=orbitals, basis=basis, normalize=normalize, eri_type = eri_type)
    elseif precision == "double"
        IntegralHelper{Float64}(molecule=molecule, orbitals=orbitals, basis=basis, normalize=normalize, eri_type = eri_type)
    else
        throw(InvalidFermiOption("precision can only be `single` or `double`. Got $precision"))
    end
end

function IntegralHelper{T}(;molecule = Molecule(), orbitals = AtomicOrbitals(), 
                           basis = Options.get("basis"), normalize = false, eri_type = nothing) where T<:AbstractFloat

    # Check if density-fitting is requested
    if Options.get("df") && eri_type === nothing
        # If the associated orbitals are AtomicOrbitals and DF is requested, JKFIT is set by default
        # Otherwise, the ERI type will be RIFIT
        eri_type = typeof(orbitals) === AtomicOrbitals ? JKFIT() : RIFIT()
    elseif !(typeof(eri_type) <: AbstractERI)
        eri_type = Chonky()
    end

    # Starts an empty cache
    cache = Dict{String,FermiMDArray{T}}()

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

# Integrals specific to orbital types
include("AtomicIntegrals.jl")
include("ROIntegrals.jl")

end # Module
