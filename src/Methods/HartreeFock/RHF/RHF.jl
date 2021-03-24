using TensorOperations
using LinearAlgebra
using Fermi.DIIS
using Fermi.Integrals: projector
using Fermi: AbstractRestrictedOrbitals

export RHF

abstract type RHFAlgorithm end

function get_scf_alg()
    implemented = [RHFa()]
    N = Options.get("scf_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(InvalidFermiOption("implementation number $N not available for RHF."))
    end
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
    C::AbstractArray{Float64,2}
end

"""
    Fermi.HartreeFock.RHF

Wave function object for Restricted Hartree-Fock methods

# High Level Interface 
```
julia> @energy rhf
# Equivalent to
julia> Fermi.HartreeFock.RHF()
```
Computes RHF using information from Fermi.CurrentOptions.

# Fields:

    molecule    Molecule object used to compute the RHF wave function
    energy      RHF Energy
    ndocc       Number of doubly occupied spatial orbitals
    nvir        Number of virtual spatial orbitals
    orbitals    RHF Orbitals

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `scf_alg`      | Picks SCF algorithm               | `String`  | [conventional]        |
| `scf_max_rms`  | RMS density convergence criterion | `Float64` | [10^-9]               |
| `scf_max_iter` | Max number of iterations          | `Int`     | [50]                  |
| `scf_e_conv`   | Energy convergence criterion      | `Float64` | [10^-10]              |
| `basis`        | What basis set to use             | `String`  | ["sto-3g"]            |
| `df`           | Whether to use density fitting    | Bool      | false                 |
| `jkfit`        | What aux. basis set to use for JK | `String`  | ["auto"]              |
| `oda`          | Whether to use ODA                | `Bool`    | [`true`]              |
| `oda_cutoff`   | When to turn ODA off (RMS)        | `Float64` | [1E-1]                |
| `oda_shutoff`  | When to turn ODA off (iter)       | `Int`     | [20]                  |
| `scf_guess`    | Which guess density to use        | `String`  | "core" ["gwh"]        |

# Lower level interfaces

    RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, ERI::Array{Float64,N}, Λ::Array{Float64,2}) where N

The RHF kernel. Computes RHF on the given `molecule` with integral information defined in `aoint`. Starts from
the given C matrix as orbitals coefficients. Λ is the orthogonalizer (S^-1/2).

_struct tree:_

**RHF** <: AbstractHFWavefunction <: AbstractWavefunction
"""
struct RHF <: AbstractHFWavefunction
    molecule::Molecule
    energy::Float64
    ndocc::Int
    nvir::Int
    orbitals::RHFOrbitals
end

function RHF(molecule::Molecule = Molecule(), ints::IntegralHelper = IntegralHelper())

    guess = Options.get("scf_guess")
    if guess == "core"
        C, Λ = RHF_core_guess(ints)
    elseif guess == "gwh"
        C, Λ = RHF_gwh_guess(molecule, ints)
    end

    RHF(molecule, ints, C, Λ, get_scf_alg())
end

function RHF(wfn::RHF)

    # Projection of A→ B done using equations described in Werner 2004 
    # https://doi.org/10.1080/0026897042000274801

    output("Using {} wave function as initial guess", wfn.orbitals.basis)

    molB = Fermi.Geometry.Molecule()

    # Assert both A and B have the same molecule.
    if molB != wfn.molecule
        throw(InvalidFermiOption(" input RHF wavefunction and Current Options have different molecules."))
    end

    basisB = Options.get("basis")
    intsB = Fermi.Integrals.IntegralHelper()

    Sbb = intsB["S"]
    Λ = Array(Sbb^(-1/2))

    Ca = wfn.orbitals.C
    Sab = projector(wfn.molecule, wfn.orbitals.basis, molB, basisB)

    T = transpose(Ca)*Sab*(Sbb^-1.0)*transpose(Sab)*Ca
    Cb = (Sbb^-1.0)*transpose(Sab)*Ca*T^(-1/2)
    Cb = real.(Cb)

    RHF(molB, intsB, FermiMDArray(Cb), FermiMDArray(Λ), get_scf_alg())
end

# Actual HF routine is in here
include("AuxRHF.jl")
# For each implementation a singleton type must be create
struct RHFa <: RHFAlgorithm end
include("RHFa.jl")
