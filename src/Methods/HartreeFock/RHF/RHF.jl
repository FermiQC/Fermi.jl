using TensorOperations
using Lints
using Fermi.DIIS
using Fermi.Integrals: ConventionalAOIntegrals, DFAOIntegrals, IntegralHelper

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

# Fields:

    molecule    Molecule object used to compute the RHF wave function
    energy      RHF Energy
    ndocc       Number of doubly occupied spatial orbitals
    nvir        Number of virtual spatial orbitals
    C           Array with MO coefficients
    eps         Array with MO energies

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
        throw(Fermi.InvalidFermiOptions("Invalid RHF algorithm: $(A)"))
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
        throw(Fermi.InvalidFermiOptions("Invalid RHF guess: $(A)"))
    end
end

"""
    Fermi.HartreeFock.RHF()

Compute RHF wave function using data from Fermi.CurrentOptions
"""
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

"""
    Fermi.HartreeFock.RHF(molecule::Molecule, aoint::IntegralHelper, )

Algorithm for to compute RHF wave function given Molecule, IntegralHelper objects.
"""
function RHF(molecule::Molecule, aoint::IntegralHelper, Alg::B, guess::GWHGuess) where B <: RHFAlgorithm 

    #form GWH guess
    @output "Using GWH Guess\n"
    S = Hermitian(aoint["S"])
    Λ = S^(-1/2)
    H = Hermitian(aoint["T"] + aoint["V"])
    ndocc = molecule.Nα#size(S,1)
    nvir = size(S,1) - ndocc
    F = Array{Float64,2}(undef, ndocc+nvir, ndocc+nvir)
    Hmax = maximum(diag(H))
    Hmin = minimum(diag(H))
    for i = 1:ndocc+nvir
        F[i,i] = H[i,i]
        for j = i+1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ*F*transpose(Λ)

    # Get orbital energies and transformed coefficients
    eps,Ct = eigen(Hermitian(Ft))

    # Reverse transformation to get MO coefficients
    C = Λ*Ct
    Co = C[:,1:ndocc]
    D = Fermi.contract(Co,Co,"um","vm")
    Eguess = RHFEnergy(D,Array(H),F)
    @output "Guess energy: {}\n" Eguess

    RHF(molecule, aoint, C, Alg)
end

function RHF(molecule::Molecule, aoint::IntegralHelper, Alg::B, guess::CoreGuess) where B <: AbstractAOIntegrals 

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

    RHF(molecule, aoint, C, Alg)
end

"""
    Fermi.HartreeFock.RHF(molecule::Molecule, aoint::IntegralHelper, Alg::B) where B <: RHFAlgorithm

Algorithm for to compute RHF wave function given Molecule, IntegralHelper objects. Starts from an input wavefunction,
and performs basis set castup (or down) as needed.
"""
function RHF(wfn::RHF, aoint::IntegralHelper, Alg::B) where B <: RHFAlgorithm 

    # Projection of A→B done using equations described in Werner 2004 
    # https://doi.org/10.1080/0026897042000274801
    @output "Using {} wave function as initial guess\n" wfn.basis
    Ca = wfn.C
    Sbb = aoint.S
    Sab = Lints.projector(wfn.LintsBasis, aoint.LintsBasis)
    T = transpose(Ca)*Sab*(Sbb^-1)*transpose(Sab)*Ca
    Cb = (Sbb^-1)*transpose(Sab)*Ca*T^(-1/2)
    Cb = real.(Cb)
    RHF(Fermi.Geometry.Molecule(), aoint, Cb, Alg)
end

function RHFEnergy(D::Array{Float64,2}, H::Array{Float64,2},F::Array{Float64,2})
    return sum(D .* (H .+ F))
end

#actual HF routine is in here
include("SCF.jl")
