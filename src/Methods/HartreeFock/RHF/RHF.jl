using TensorOperations
using Lints
using Fermi.DIIS
using Fermi.Integrals: IntegralHelper
using Fermi.Orbitals: CanonicalOrbital, CanonicalOrbitals, OrbDict

# Define Algorithims
abstract type RHFAlgorithm end
struct ConventionalRHF <: RHFAlgorithm end
struct DFRHF <: RHFAlgorithm end

# Define Guesses
abstract type RHFGuess end
struct CoreGuess   <: RHFGuess end
struct GWHGuess    <: RHFGuess end
struct HuckelGuess <: RHFGuess end

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
| `scf_alg`      | picks SCF algorithm               | `String`  | "df" ["conventional"] |
| `scf_max_rms`  | RMS density convergence criterion |`Float64`  | [10^-10]              |
| `scf_max_iter` | Max number of iterations          | `Int`     | [50]                  |
| `basis`        | What basis set to use             | `String`  | ["sto-3g"]            |
| `jkfit`        | What aux. basis set to use for JK | `String`  | ["auto"]              |
| `oda`          | Whether to use ODA                | `Bool`    | [`True`]              |
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
    eps::Array{Float64,1}
    ints::IntegralHelper
end


function select_alg(A::String)
    implemented = Dict{String,Any}(
        "conventional" => (ConventionalRHF()),
        "df"           => (DFRHF())
       )

    try
        return implemented[A]
    catch KeyError
        throw(Fermi.InvalidFermiOption("Invalid RHF algorithm: $(A)"))
    end
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
    Alg = select_alg(Fermi.CurrentOptions["scf_alg"])
    ints = Fermi.Integrals.IntegralHelper()
    guess = select_guess(Fermi.CurrentOptions["scf_guess"])
    RHF(molecule, ints, Alg, guess)
end

function RHF(molecule::Molecule, ints::IntegralHelper, Alg::B, guess::GWHGuess) where B <: RHFAlgorithm 

    #form GWH guess
    @output "Using GWH Guess\n"
    S = ints["S"]
    d, U = diagonalize(S, sortby = x->1/abs(x))
    Λ = S^(-1/2)
    idxs = [abs(d[i]) > 1E-7 for i=1:size(S,1)]
    @output "Found {} linear dependencies. Projected them out.\n" size(S,1) - sum(idxs)
    #Λ = FermiMDArray(real.(Λ[:,idxs]))
    H = Hermitian(ints["T"] + ints["V"])
    ndocc = molecule.Nα
    nvir = size(S,1) - ndocc
    F = similar(S)
    Hmax = maximum(diag(H))
    Hmin = minimum(diag(H))
    for i = 1:ndocc+nvir
        F[i,i] = H[i,i]
        for j = i+1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ'*F*Λ

    # Get orbital energies and transformed coefficients
    eps, Ct = diagonalize(Hermitian(Ft))

    # Reverse transformation to get MO coefficients
    C = Λ*Ct
    Co = C[:,1:ndocc]
    ######
    #D = Fermi.contract(Co,Co,"um","vm")
    D = @tensor D[u,v] := Co[u,m]*Co[v,m]
    Eguess = RHFEnergy(D, H, F)
    @output "Guess energy: {}\n" Eguess

    RHF(molecule, ints, C, Λ, Alg)
end

function RHF(molecule::Molecule, aoint::IntegralHelper, Alg::B, guess::CoreGuess) where B <: RHFAlgorithm

    #Form core guess
    @output "Using Core Guess\n"
    S = Hermitian(aoint["S"])
    Λ = S^(-1/2)
    H = Hermitian(aoint["T"] + aoint["V"])
    F = Array{Float64,2}(undef, ndocc+nvir, ndocc+nvir)
    F .= H
    Ft = Λ*F*transpose(Λ)

    # Get orbital energies and transformed coefficients
    eps,Ct = eigen(Hermitian(Ft))

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    RHF(molecule, aoint, C, Λ, Alg)
end

function RHF(wfn::RHF)
    Alg = select_alg(Fermi.CurrentOptions["scf_alg"])
    ints = Fermi.Integrals.IntegralHelper()
    RHF(wfn, ints, Alg)
end

function RHF(wfn::RHF, aoint::IntegralHelper, Alg::B) where B <: RHFAlgorithm 

    # Projection of A→ B done using equations described in Werner 2004 
    # https://doi.org/10.1080/0026897042000274801
    @output "Using {} wave function as initial guess\n" wfn.ints.bname["primary"]

    Ca = Float64[]
    for orb in wfn.ints.orbs["FU"].orbs
        Ca = vcat(Ca, orb.C)
    end

    nbf = Int(√length(Ca))
    Ca = reshape(Ca, (nbf, nbf))

    Sbb = aoint["S"]
    S = Hermitian(aoint["S"])
    Λ = Array(S^(-1/2))

    open("/tmp/molfile1.xyz","w") do molfile
        natom = length(wfn.molecule.atoms)
        write(molfile,"$natom\n\n")
        write(molfile,Fermi.Geometry.get_xyz(wfn.molecule))
    end

    Lints.libint2_init()
    wfnmol = Lints.Molecule("/tmp/molfile1.xyz")
    wfnbas = Lints.BasisSet(wfn.ints.bname["primary"], wfnmol)

    molecule = Fermi.Geometry.Molecule()
    open("/tmp/molfile2.xyz","w") do molfile
        natom = length(molecule.atoms)
        write(molfile,"$natom\n\n")
        write(molfile,Fermi.Geometry.get_xyz(molecule))
    end

    mol = Lints.Molecule("/tmp/molfile2.xyz")
    bas = Lints.BasisSet(aoint.bname["primary"], wfnmol)

    Sab = Lints.projector(wfnbas, bas)
    T = transpose(Ca)*Sab*(Sbb^-1)*transpose(Sab)*Ca
    Cb = (Sbb^-1)*transpose(Sab)*Ca*T^(-1/2)
    Cb = real.(Cb)
    RHF(molecule, aoint, Cb, Λ, Alg)
end

function RHFEnergy(D::FermiMDArray{Float64}, H::FermiMDArray{Float64},F::FermiMDArray{Float64})
    return sum(D .* (H .+ F))
end

#actual HF routine is in here
include("SCF.jl")
