using Fermi.DIIS

using TensorOperations
using LinearAlgebra
using Formatting
import Base: show

export UHF

abstract type UHFAlgorithm end

struct UHF <: AbstractHFWavefunction
    molecule::Molecule
    energy::Float64
    ndocc::Int
    nsocc::Int
    nvir::Int
    orbitals::UHFOrbitals
    e_conv::Float64
    d_conv::Float64
end

function UHF(x...)
    if !any(i-> i isa UHFAlgorithm, x)
        UHF(x..., get_scf_alg())
    else
        throw(FermiException("invalid arguments for UHF method: $(x[1:end-1])"))
    end
end

struct UHFa <: UHFAlgorithm end
include("UHFa.jl")

function string_repr(X::UHF)
    out = ""
    out = out*" ⇒ Fermi Unrestricted Hartree--Fock Wave function\n"
    out = out*" ⋅ Basis:                  $(X.orbitals.basis)\n"
    out = out*" ⋅ Energy:                 $(X.energy)\n"
    out = out*" ⋅ Occ. Spatial Orbitals:  $(X.ndocc)\n"
    out = out*" ⋅ Vir. Spatial Orbitals:  $(X.nvir)\n"
    out = out*"Convergence: " 
    out = out*"ΔE => $(format("{:1.2e}",abs(X.e_conv)))"
    out = out*" Dᵣₘₛ => $(format("{:1.2e}",abs(X.d_conv)))"
    return out
end

function show(io::IO, ::MIME"text/plain", X::UHF)
    print(string_repr(X))
end


