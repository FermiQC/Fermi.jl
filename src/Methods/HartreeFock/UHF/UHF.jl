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
    nocc::Int
    nvir::Int
    orbitals::UHFOrbitals
    e_conv::Float64
    d_conv::Float64
end

function UHF(x...)
    if !any(i-> i isa UHFAlgorithm, x)
        UHF(x..., get_uhf_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for UHF method: $args"))
    end
end

"""
    Fermi.HartreeFock.get_uhf_alg()

Returns a singleton type corresponding to a UHF implementation.
"""
function get_uhf_alg(N::Int = Options.get("uhf_alg"))
    try 
        return get_uhf_alg(Val(N))
    catch MethodError
        throw(FermiException("implementation number $N not available for UHF."))
    end
end

struct UHFa <: UHFAlgorithm end
include("UHFa.jl")
include("UHFHelper.jl")
# And a number is assigned to the implementation
get_uhf_alg(x::Val{1}) = UHFa()

function string_repr(X::UHF)
    out = ""
    out = out*" ⇒ Fermi Unrestricted Hartree--Fock Wave function\n"
    out = out*" ⋅ Basis:                  $(X.orbitals.basis)\n"
    out = out*" ⋅ Energy:                 $(X.energy)\n"
    out = out*" ⋅ Occ. Spatial Orbitals:  $(X.nocc)\n"
    out = out*" ⋅ Vir. Spatial Orbitals:  $(X.nvir)\n"
    out = out*"Convergence: " 
    out = out*"ΔE => $(format("{:1.2e}",abs(X.e_conv)))"
    out = out*" Dᵣₘₛ => $(format("{:1.2e}",abs(X.d_conv)))"
    return out
end

function show(io::IO, ::MIME"text/plain", X::UHF)
    print(string_repr(X))
end


