using TensorOperations
using Lints
using Fermi.DIIS
using Fermi.Integrals: ConventionalAOIntegrals, DFAOIntegrals

# Define Algorithims
abstract type RHFAlgorithm end
struct ConventionalRHF <: RHFAlgorithm end
struct DFRHF <: RHFAlgorithm end

# Define Guesses
abstract type RHFGuess end
struct CoreGuess <: RHFGuess end
struct GWHGuess  <: RHFGuess end

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
    basis::String
#    LintsBasis::Lints.BasisSetAllocated
    molecule::Molecule
    energy::Float64
    ndocc::Int
    nvir::Int
    C::Array{Float64,2} 
    eps::Array{Float64,1}
    ints::I where I <: AbstractAOIntegrals
end

# Algorithm-specific dispatches
#include("ConventionalRHF.jl")
#include("DF-RHF.jl")
include("SCF.jl")

function select_alg(A::String)
    implemented = Dict{String,Any}(
        "conventional" => (ConventionalAOIntegrals,ConventionalRHF()),
        "df"           => (DFAOIntegrals,DFRHF())
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

# General dispatches
"""
    Fermi.HartreeFock.RHF()

Compute RHF wave function using data from Fermi.CurrentOptions
"""
function RHF()
    molecule = Molecule()
    RHF(molecule)
end

function RHF(molecule::Molecule)
    aoint_type,Alg = select_alg(Fermi.CurrentOptions["scf_alg"])
    aoint = aoint_type(molecule)
    guess = select_guess(Fermi.CurrentOptions["scf_guess"])
    RHF(molecule, aoint, Alg, guess)
end

"""
    Fermi.HartreeFock.RHF(molecule::Molecule, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

Conventional algorithm for to compute RHF wave function given Molecule, Integrals objects.
"""
function RHF(molecule::Molecule, aoint::A, Alg::B, guess::GWHGuess) where { A <: AbstractAOIntegrals,
                                                                            B <: RHFAlgorithm }

    @output "Using GWH Guess\n"
    S = Hermitian(aoint.S)
    Λ = S^(-1/2)
    H = Hermitian(aoint.T + aoint.V)
    ndocc = molecule.Nα#size(S,1)
    nvir = size(S,1) - ndocc
    F = Array{Float64,2}(undef, ndocc+nvir, ndocc+nvir)
    for i = 1:ndocc+nvir
        F[i,i] = H[i,i]
        for j = 1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ*F*transpose(Λ)

    # Get orbital energies and transformed coefficients
    eps,Ct = eigen(Hermitian(Ft))

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    RHF(molecule, aoint, C, Alg)
end

function RHF(molecule::Molecule, aoint::A, Alg::B, guess::CoreGuess) where { A <: AbstractAOIntegrals,
                                                                             B <: AbstractAOIntegrals }

    S = Hermitian(aoint.S)
    Λ = S^(-1/2)
    H = Hermitian(aoint.T + aoint.V)
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
    Fermi.HartreeFock.RHF(wfn::RHF, aoint::ConventionalAOIntegrals, Alg::ConventionalRHF)

Conventional algorithm for to compute RHF wave function. Inital guess for orbitals is built from given RHF wfn. Integrals
are taken from the aoint input.
"""
function RHF(wfn::RHF, aoint::A, Alg::B) where { A <: AbstractAOIntegrals,
                                                 B <: RHFAlgorithm }

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

"""
    Fermi.HartreeFock.RHF(wfn::RHF)

    Compute RHF wave function using the input RHF wave function (wfn) to generate a guess for orbitals.
"""
function RHF(wfn::RHF)
    aoint_type,Alg = select_algorithm(Fermi.CurrentOptions["scf_alg"])
    RHF(wfn, Alg)
end

function RHFEnergy(D::Array{Float64,2}, H::Array{Float64,2},F::Array{Float64,2})
    return sum(D .* (H .+ F))
end
