using TensorOperations
using LinearAlgebra
using Lints
using Fermi.DIIS
using Fermi.Integrals: IntegralHelper, projector
using Fermi: AbstractOrbitals, Options
using Fermi.Geometry: Molecule

# Define Guesses
abstract type RHFGuess end
struct CoreGuess   <: RHFGuess end
struct GWHGuess    <: RHFGuess end

# Define RHF Orbital

struct RHFOrbitals <: AbstractOrbitals
    molecule::Molecule
    basis::String
    aux::String
    eps::AbstractArray{Float64,1}
    C::AbstractArray{Float64,2}
end

"""
    Fermi.HartreeFock.RHF

Wave function object for Restricted Hartree-Fock methods

# High Level Interface 
    RHF()

Computes RHF using information from Fermi.CurrentOptions.

# Fields:

    molecule    Molecule object used to compute the RHF wave function
    energy      RHF Energy
    ndocc       Number of doubly occupied spatial orbitals
    nvir        Number of virtual spatial orbitals
    C           Array with MO coefficients
    eps         Array with MO energies

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `scf_type`     | picks SCF algorithm               | `String`  | "df" ["conventional"] |
| `scf_max_rms`  | RMS density convergence criterion |`Float64`  | [10^-10]              |
| `scf_max_iter` | Max number of iterations          | `Int`     | [50]                  |
| `basis`        | What basis set to use             | `String`  | ["sto-3g"]            |
| `jkfit`        | What aux. basis set to use for JK | `String`  | ["auto"]              |
| `oda`          | Whether to use ODA                | `Bool`    | [`true`]              |
| `oda_cutoff`   | When to turn ODA off (RMS)        | `Float64` | [1E-1]                |
| `oda_shutoff`  | When to turn ODA off (iter)       | `Int`     | [20]                  |
| `scf_guess`    | Which guess density to use        |           | "core" ["gwh"]        |

# Lower level interfaces

    RHF(molecule::Molecule, aoint::IntegralHelper, C::Array{Float64,2}, ERI::Array{Float64,N}, Λ::Array{Float64,2}) where N

The RHF kernel. Computes RHF on the given `molecule` with integral information defined in `aoint`. Starts from
the given C matrix. 

_struct tree:_

**RHF** <: AbstractHFWavefunction <: AbstractReferenceWavefunction <: AbstractWavefunction
"""
struct RHF <: AbstractHFWavefunction
    molecule::Molecule
    energy::Float64
    ndocc::Int
    nvir::Int
    orbitals::RHFOrbitals
end

function select_guess(A::String)
    implemented = Dict{String,Any}(
        "gwh"  => GWHGuess(),
        "core" => CoreGuess()
    )
    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOption("Invalid RHF guess: $(A)"))
    end
end

function RHF()
    molecule = Molecule()
    RHF(molecule)
end

function RHF(molecule::Molecule)
    guess = select_guess(Options.get("scf_guess"))
    ints = Fermi.Integrals.IntegralHelper()
    RHF(molecule, ints, guess)
end

function RHF(molecule::Molecule, ints::IntegralHelper, guess::GWHGuess)

    # Form GWH guess
    output("Using GWH Guess")
    S = ints["S"]
    d, U = diagonalize(S, sortby = x->1/abs(x))
    Λ = S^(-1/2)
    idxs = [abs(d[i]) > 1E-7 for i = eachindex(d)]

    output("Found {} linear dependencies", length(d) - sum(idxs))

    H = real.(ints["T"] + ints["V"])
    ndocc = molecule.Nα
    nvir = size(S,1) - ndocc
    F = similar(S)

    for i = 1:ndocc+nvir
        F[i,i] = H[i,i]
        for j = i+1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ'*F*Λ

    # Get orbital energies and transformed coefficients
    eps, Ct = diagonalize(Ft, hermitian=true)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct
    Co = C[:,1:ndocc]

    @tensor D[u,v] := Co[u,m]*Co[v,m]

    Eguess = RHFEnergy(D, H, F)

    RHF(molecule, ints, C, Λ)
end

function RHF(molecule::Molecule, ints::IntegralHelper, guess::CoreGuess)

    # Form Core Guess
    output("Using Core Guess")
    S = ints["S"]
    Λ = S^(-1/2)
    F = ints["T"] + ints["V"]
    Ft = Λ*F*Λ'

    # Get orbital energies and transformed coefficients
    _, Ct = diagonalize(Ft, hermitian=true)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    RHF(molecule, ints, C, Λ, Alg)
end

function RHF(wfn::RHF)

    # Start RHF computation from another RHF object
    ints = Fermi.Integrals.IntegralHelper(wfn.molecule, wfn.orbitals.basis, wfn.orbitals.aux)
    RHF(wfn, ints)
end

function RHF(wfn::RHF, Bints::IntegralHelper)

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

    RHF(molB, intsB, FermiMDArray(Cb), FermiMDArray(Λ))
end

include("AuxRHF.jl")
# Actual HF routine is in here
include("SCF.jl")
