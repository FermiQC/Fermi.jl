using Fermi.HartreeFock

export RCCSD

abstract type RCCSDAlgorithm end

function get_rccsd_alg()
    implemented = [RCCSDa()]
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
    energy::T
    correlation::T
    T1::AbstractArray{T,2}
    T2::AbstractArray{T,4}
    δE::T
    residue::T
end

function RCCSD()
    aoints = IntegralHelper{Float64}()
    rhf = RHF(aoints)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RCCSD(moints, aoints)
end

function RCCSD(moints::IntegralHelper{T1,Chonky,O}, aoints::IntegralHelper{T2,Chonky,AtomicOrbitals}) where {T1<:AbstractFloat,
                                                                                        T2<:AbstractFloat,O<:AbstractOrbitals}
    Fermi.Integrals.mo_from_ao(moints, aoints, "Fd","OOOO", "OOOV", "OOVV", "OVOV", "OVVV", "VVVV")
    RCCSD(moints)
end

function RCCSD(moints::IntegralHelper{T1,E1,O}, aoints::IntegralHelper{T2,E2,AtomicOrbitals}) where {T1<:AbstractFloat,T2<:AbstractFloat,
                                                                                E1<:AbstractDFERI,E2<:AbstractDFERI,O<:AbstractOrbitals}
    Fermi.Integrals.mo_from_ao(moints, aoints, "Fd","BOO", "BOV", "BVV")
    RCCSD(moints)
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

function RCCSD(moints::IntegralHelper{T,E,O}) where {T<:AbstractFloat,E<:AbstractERI,O<:AbstractRestrictedOrbitals}

    # Create zeroed guesses for amplitudes

    o = moints.molecule.Nα - Options.get("drop_occ")
    v = size(moints.orbitals.C,1) - Options.get("drop_vir") - moints.molecule.Nα

    output("Using MP2 guess")
    T1guess = moints["Fia"]
    T2guess = permutedims(moints["OVOV"], (1,3,2,4))

    # Orbital energies line
    if haskey(moints.cache, "D1")
        d = moints["D1"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        frozen = Options.get("drop_occ")
        inac = Options.get("drop_vir")
        ϵo = Fd[(1+frozen:ndocc)]
        ϵv = Fd[(1+ndocc):end-inac]

        d = FermiMDArray([ϵo[i]-ϵv[a] for i=eachindex(ϵo), a=eachindex(ϵv)])
        moints["D1"] = d
    end

    if haskey(moints.cache, "D2")
        D = moints["D2"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        frozen = Options.get("drop_occ")
        inac = Options.get("drop_vir")
        ϵo = Fd[(1+frozen:ndocc)]
        ϵv = Fd[(1+ndocc):end-inac]

        D = FermiMDArray([ϵo[i]+ϵo[j]-ϵv[a]-ϵv[b] for i=eachindex(ϵo), j=eachindex(ϵo), a=eachindex(ϵv), b=eachindex(ϵv)])
        moints["D2"] = D
    end

    T1guess ./= d
    T2guess ./= D

    RCCSD(moints, T1guess, T2guess, get_rccsd_alg())
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
struct RCCSDa <: RCCSDAlgorithm end
include("RCCSDa.jl")