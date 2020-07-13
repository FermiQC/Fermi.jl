"""
    Fermi.Integrals

Module to compute integrals using Lints.jl
"""
module Integrals

using Fermi
using Fermi.Geometry: Molecule
using Lints
using LinearAlgebra

export AbstractIntegrals
export AOIntegrals

abstract type AbstractIntegrals end
abstract type AbstractAOIntegrals  <: AbstractIntegrals end
abstract type AbstractMOIntegrals  <: AbstractIntegrals end
abstract type AbstractRMOIntegrals <: AbstractIntegrals end

struct ConventionalAOIntegrals{T} <: AbstractAOIntegrals where T <: AbstractFloat
    S::Array{T,2}
    T::Array{T,2}
    V::Array{T,2}
    ERI:: I where I <: Fermi.AbstractTensor
end

function ConventionalAOIntegrals()
    basis = Fermi.CurrentOptions["basis"]
    ConventionalAOIntegrals(Molecule(), basis)
end

function ConventionalAOIntegrals(molecule::Molecule)
    basis = Fermi.CurrentOptions["basis"]
    ConventionalAOIntegrals(molecule, basis)
end

function ConventionalAOIntegrals(basis::String)

    ConventionalAOIntegrals(Molecule(), basis)
end

function ConventionalAOIntegrals(molecule::Molecule, basis::String)

    ConventionalAOIntegrals(molecule, basis, Fermi.ComputeEnvironment.interconnect,
                                 Fermi.ComputeEnvironment.communicator,
                                 Fermi.ComputeEnvironment.accelerator)
end

function ConventionalAOIntegrals(molecule::Molecule, basis::String, interconnect::Fermi.Environments.No_IC,
                                                        communicator::Fermi.Environments.NoCommunicator,
                                                        accelerator::Fermi.Environments.NoAccelerator)

    open("/tmp/molfile.xyz","w") do molfile
        natom = length(molecule.atoms)
        write(molfile,"$natom\n\n")
        write(molfile,Fermi.Geometry.get_xyz(molecule))
    end

    Lints.libint2_init()
    mol = Lints.Molecule("/tmp/molfile.xyz")
    bas = Lints.BasisSet(basis, mol)

    nprim = Lints.max_nprim(bas)
    l = Lints.max_l(bas)
    S_engine = Lints.OverlapEngine(nprim,l)
    T_engine = Lints.KineticEngine(nprim,l)
    V_engine = Lints.NuclearEngine(nprim,l,mol)
    I_engines = []
    sz = Lints.getsize(S_engine,bas)
    for i in 1:Threads.nthreads()
        push!(I_engines,Lints.ERIEngine(nprim,l))
    end
    S = zeros(sz,sz)
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    I = zeros(sz,sz,sz,sz)
    Lints.make_2D(S,S_engine,bas)
    Lints.make_2D(T,T_engine,bas)
    Lints.make_2D(V,V_engine,bas)
    Lints.make_ERI(I,I_engines,bas)
    I = Fermi.MemTensor(I)
    Lints.libint2_finalize()

    return ConventionalAOIntegrals{Float64}(S, T, V, I)
end

struct PhysMOIntegrals{T} <: AbstractRMOIntegrals where T <: AbstractFloat
   oooo::Array{T,4}
   ooov::Array{T,4}
   oovv::Array{T,4}
   ovov::Array{T,4}
   ovvv::Array{T,4}
   vvvv::Array{T,4}
   oo::Array{T,2}
   ov::Array{T,2}
   vv::Array{T,2}
end

end #module
