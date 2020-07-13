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

"""
    Fermi.Integrals.AbstractIntegrals

Abstract type common to all integral objects.

_struct tree:_

**AbstractIntegrals** (Top level)
"""
abstract type AbstractIntegrals end
"""
    Fermi.Integrals.AbstractAOIntegrals

Abstract type common to all *A*tomic *O*rbitals integral objects.

_struct tree:_

**AbstractAOIntegrals** <: AbstractIntegrals
"""
abstract type AbstractAOIntegrals  <: AbstractIntegrals end
"""
    Fermi.Integrals.AbstractMOIntegrals

Abstract type common to all *M*tomic *O*rbitals integral objects.

_struct tree:_

**AbstractMOIntegrals** <: AbstractIntegrals
"""
abstract type AbstractMOIntegrals  <: AbstractIntegrals end

"""
    Fermi.Integrals.ConventionalAOIntegrals

Object holding AO integrals in memory.

_struct tree:_

**ConventionalAOIntegrals** <: AbstractAOIntegrals <: AbstractIntegrals
"""
struct ConventionalAOIntegrals{T} <: AbstractAOIntegrals where T <: AbstractFloat
    S::Array{T,2}
    T::Array{T,2}
    V::Array{T,2}
    ERI:: I where I <: Fermi.AbstractTensor
end

"""
    Fermi.Integrals.ConventionalAOIntegrals()

Uses data from Fermi.CurrentOptions to compute and return a ConventionalAOIntegrals object.
"""
function ConventionalAOIntegrals()
    basis = Fermi.CurrentOptions["basis"]
    ConventionalAOIntegrals(Molecule(), basis)
end

"""
    Fermi.Integrals.ConventionalAOIntegrals(molecule::Molecule)

Uses the molecule object and basis from Fermi.CurrentOptions to compute and return a 
ConventionalAOIntegrals object.
"""
function ConventionalAOIntegrals(molecule::Molecule)
    basis = Fermi.CurrentOptions["basis"]
    ConventionalAOIntegrals(molecule, basis)
end

"""
    Fermi.Integrals.ConventionalAOIntegrals(basis::String)

For the given basis, use the Molecule store in Fermi.CurrentOptions to compute and return
a ConventionalAOIntegrals object.
"""
function ConventionalAOIntegrals(basis::String)

    ConventionalAOIntegrals(Molecule(), basis)
end

"""
    Fermi.Integrals.ConventionalAOIntegrals(molecule::Molecule, basis::String)

For the given molecule and basis, uses the computing environment from Fermi.ComputeEnvironment
to compute and return a ConventionalAOIntegrals object.
"""
function ConventionalAOIntegrals(molecule::Molecule, basis::String)

    ConventionalAOIntegrals(molecule, basis, Fermi.ComputeEnvironment.interconnect,
                                 Fermi.ComputeEnvironment.communicator,
                                 Fermi.ComputeEnvironment.accelerator)
end

"""
    Fermi.Integrals.ConventionalAOIntegrals(molecule::Molecule, basis::String, 
                                                  interconnect::Fermi.Environments.No_IC,
                                                  communicator::Fermi.Environments.NoCommunicator,
                                                  accelerator::Fermi.Environments.NoAccelerator)

For the given molecule, basis, and comuting environment compute and return a
ConventionalAOIntegrals object.
"""
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

"""
    Fermi.Integrals.PhysRMOIntegrals
Object holding restricted MO integrals in memory using Physicists' notation.

_struct tree:_

**PhysRMOIntegrals** <: AbstractMOIntegrals <: AbstractIntegrals
"""
struct PhysRMOIntegrals{T} <: AbstractMOIntegrals where T <: AbstractFloat
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
