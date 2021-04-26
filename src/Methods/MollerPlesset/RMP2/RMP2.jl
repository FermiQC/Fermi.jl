using Fermi.HartreeFock: RHF, RHFOrbitals

export RMP2

abstract type RMP2Algorithm end

function get_mp2_alg()
    implemented = [RMP2a()]
    N = Options.get("mp2_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(InvalidFermiOption("implementation number $N not available for RMP2."))
    end
end

"""
    RMP2{T} <: AbstractMPWavefunction

    TODO
"""
struct RMP2{T} <: AbstractMPWavefunction
    correlation::T
    energy::T
end

function RMP2()
    aoints = IntegralHelper{Float64}()
    rhf = RHF(aoints)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RMP2(moints, aoints)
end

function RMP2(O::AbstractRestrictedOrbitals)
    aoints = IntegralHelper()
    moints = IntegralHelper(orbitals=O)
    RMP2(moints, aoints)
end

function RMP2(rhf::RHF)
    aoints = IntegralHelper()
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RMP2(moints, aoints)
end

function RMP2(M::Molecule)
    aoints = IntegralHelper{Float64}(molecule=M)
    rhf = RHF(M,aoints)
    moints = IntegralHelper(molecule = M, orbitals=rhf.orbitals)
    RMP2(moints, aoints)
end

function RMP2(moints::IntegralHelper{T1,Chonky,O}, aoints::IntegralHelper{T2,Chonky,AtomicOrbitals}) where {T1<:AbstractFloat,
                                                                                        T2<:AbstractFloat,O<:AbstractOrbitals}
    Fermi.Integrals.mo_from_ao(moints, aoints, "Fia","OVOV")
    RMP2(moints)
end

function RMP2(moints::IntegralHelper{T1,E1,O}, aoints::IntegralHelper{T2,E2,AtomicOrbitals}) where {T1<:AbstractFloat,T2<:AbstractFloat,
                                                                                E1<:AbstractDFERI,E2<:AbstractDFERI,O<:AbstractOrbitals}
    Fermi.Integrals.mo_from_ao(moints, aoints, "Fia","BOV")
    RMP2(moints)
end

function RMP2(ints::IntegralHelper{T,E,O}) where {T<:AbstractFloat,E<:AbstractERI,O<:AbstractRestrictedOrbitals}
    RMP2(ints, get_mp2_alg())
end

# For each implementation a singleton type must be create
struct RMP2a <: RMP2Algorithm end
include("RMP2a.jl")