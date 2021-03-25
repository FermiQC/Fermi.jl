using Fermi.HartreeFock.RHF
using Fermi: AbstractRestrictedOrbitals

export RCCSD

abstract type RCCSDAlgorithm end

function get_rccsd_alg()
    implemented = [RCCSD1()]
    N = Options.get("cc_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(InvalidFermiOption("implementation number $N not available for RCCSD."))
    end
end

"""
    Fermi.CoupledCluster.RCCSD
    TODO

Fermi struct that holds information about RCCSD wavefunctions

# Fields

    CorrelationEnergy   CCSD correlation energy
    T1                  T1 amplitudes
    T2                  T2 amplitudes

_struct tree:_

**RCCSD** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
struct RCCSD{T} <: AbstractCCWavefunction 
    guessenergy::T
    correlation::T
    T1::AbstractArray{T,2}
    T2::AbstractArray{T,4}
    converged::Bool
end

function RCCSD(x...)
    precision = Options.get("precision")
    if precision == "single"
        RCCSD{Float32}(x...)
    elseif precision == "double"
        RCCSD{Float64}(x...)
    else
        throw(InvalidFermiOption("precision can only be `single` or `double`. Got $precision"))
    end
end

function RCCSD{Float64}(mol = Molecule(), ints = IntegralHelper{Float64}())

    # Compute Restricted Hartree-Fock
    refwfn = RHF(mol, ints)

    # Delete integrals that are not gonna be used anymore
    delete!(ints, "S", "T", "V", "JKERI")
    RCCSD{Float64}(refwfn, ints)
end

function RCCSD{Float32}(mol = Molecule(), ints = IntegralHelper{Float32}())

    # Note that using this method a new IntegralHelper object
    # is created within the RHF call. This is necessary becasue
    # RHF needs a double precision integral helper

    # Compute Restricted Hartree-Fock
    refwfn = RHF(mol)

    RCCSD{Float32}(refwfn, ints)
end

function RCCSD{T}(refwfn::RHF, ints::IntegralHelper{T})

    # Create zeroed guesses for amplitudes

    o = refwfn.ndocc - Options.get("drop_occ")
    v = refwfn.nvir - Options.get("drop_nvir")
    T1guess = FermiMDzeros(T, o, v)
    T2guess = FermiMDzeros(T, o, o, v, v)
    RCCSD{T}(refwfn, ints, T1guess, T2guess, alg = get_rccsd_alg())
end


"""
    Fermi.CoupledCluster.RCCSD{T}()

Compute a RCCSD wave function for a given precision T (Float64 or Float32)
"""
function RCCSD{T}(guess::RCCSD{Tb}) where { T <: AbstractFloat,
                                           Tb <: AbstractFloat }
    alg = select_algorithm(Fermi.CurrentOptions["cc_alg"])
    RCCSD{T}(guess,alg)
end

# For each implementation a singleton type must be create
struct RCCSD1 <: RCCSDAlgorithm end
include("RCCSDa.jl")