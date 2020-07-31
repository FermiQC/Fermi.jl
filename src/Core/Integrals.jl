"""
    Fermi.Integrals

Module to compute integrals using Lints.jl
"""
module Integrals

using Fermi
using Fermi.Output
using Fermi.Geometry: Molecule
using Lints
using LinearAlgebra
using TensorOperations

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

# Fields:

    S    Overlap AO matrix
    T    Kinetic energy AO matrix
    V    Nuclear attraction matrix
    ERI  Tensor with electron repulsion integrals

_struct tree:_

**ConventionalAOIntegrals** <: AbstractAOIntegrals <: AbstractIntegrals
"""
struct ConventionalAOIntegrals{T} <: AbstractAOIntegrals where T <: AbstractFloat
    basis::String
    LintsBasis::Lints.BasisSetAllocated
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
    @output "   • Lints started\n\n"
    @output "   Basis set: {}\n" basis
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
    #Lints.libint2_finalize()
    @output "Exiting Lints.\n\n"

    return ConventionalAOIntegrals{Float64}(basis, bas, S, T, V, I)
end

struct DFAOIntegrals{T} <: AbstractAOIntegrals where T <: AbstractFloat
    bname::String
    dfbname::String
    basis::Lints.BasisSetAllocated
    dfbasis::Lints.BasisSetAllocated
    S::Array{T,2}
    T::Array{T,2}
    V::Array{T,2}
    B::Array{T,3}
end

function DFAOIntegrals(mol::Fermi.Geometry.Molecule)
    bname = Fermi.CurrentOptions["basis"]
    DFAOIntegrals(mol,bname)
end
function DFAOIntegrals()
    bname = Fermi.CurrentOptions["basis"]
    DFAOIntegrals(Molecule(),bname)
end

function DFAOIntegrals(molecule::Molecule, bname::String)
    DFAOIntegrals(molecule,bname,Fermi.ComputeEnvironment.interconnect,
                            Fermi.ComputeEnvironment.communicator,
                            Fermi.ComputeEnvironment.accelerator)
end

function DFAOIntegrals(molecule::Molecule, bname::String, interconnnect::Fermi.Environments.No_IC,
                                                communicator::Fermi.Environments.NoCommunicator,
                                                accelerator::Fermi.Environments.NoAccelerator)
    open("/tmp/molfile.xyz","w") do molfile
        natom = length(molecule.atoms)
        write(molfile,"$natom\n\n")
        write(molfile,Fermi.Geometry.get_xyz(molecule))
    end

    @assert lowercase(bname) in ["cc-pvdz",
                                  "cc-pvtz",
                                  "cc-pvqz",
                                  "cc-pv5z",
                                  "aug-cc-pvdz",
                                  "aug-cc-pvtz",
                                  "aug-cc-pvqz",
                                  "aug-cc-pv5z"] "Only Dunning basis sets are supported for density fitting!"
    dfbname = bname*"-RI"
    Lints.libint2_init()
    @output "   • Lints started\n\n"
    @output "   Basis set: {}\n" bname
    @output "   DF Basis set: {}\n" dfbname 
    mol = Lints.Molecule("/tmp/molfile.xyz")
    bas = Lints.BasisSet(bname, mol)
    dfbas = Lints.BasisSet(dfbname,mol)

    nprim = max(Lints.max_nprim(bas),Lints.max_nprim(dfbas))
    l = max(Lints.max_l(bas),Lints.max_l(dfbas))

    S_engine = Lints.OverlapEngine(nprim,l)
    T_engine = Lints.KineticEngine(nprim,l)
    V_engine = Lints.NuclearEngine(nprim,l,mol)
    eri_engines = [Lints.DFEngine(nprim,l) for i=1:Threads.nthreads()]
    sz = Lints.getsize(S_engine,bas)
    dfsz = Lints.getsize(S_engine,dfbas)
    S = zeros(sz,sz)
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    J = zeros(dfsz,dfsz)
    Pqp = zeros(dfsz,sz,sz)
    Lints.make_2D(S,S_engine,bas)
    Lints.make_2D(T,T_engine,bas)
    Lints.make_2D(V,V_engine,bas)
    Lints.make_j(J,eri_engines[1],dfbas)
    Lints.make_b(Pqp,eri_engines,bas,dfbas)
    Jh = J^(-1/2)
    B = zeros(dfsz,sz,sz)
    Fermi.contract!(B,Pqp,Jh,"Qpq","Pqp","PQ")
    @output "Exiting Lints.\n\n"

    return DFAOIntegrals{Float64}(bname,dfbname,bas, dfbas, S, T, V, B)
end

"""
    Fermi.Integrals.PhysRMOIntegrals
Object holding restricted MO integrals in memory using Physicists' notation.

_struct tree:_

**PhysRMOIntegrals** <: AbstractMOIntegrals <: AbstractIntegrals
"""
struct PhysRestrictedMOIntegrals{T} <: AbstractMOIntegrals where T <: AbstractFloat
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

"""
test
"""
function PhysRestrictedMOIntegrals{T}(ndocc::Int, nvir::Int, drop_occ::Int, drop_vir::Int, C::Array{Float64,2}, aoint::ConventionalAOIntegrals) where T <: AbstractFloat

    nmo = ndocc + nvir
    Co = C[:, (1+drop_occ):ndocc]
    Cv = C[:, (1+ndocc):(nmo-drop_vir)]

    # Get MO ERIs

    TensorType = typeof(aoint.ERI)
    ERI = TensorType(permutedims(aoint.ERI.data,(1,3,2,4)))
    oooo = transform_eri(ERI, Co, Co, Co, Co)
    oooo = T.(oooo)

    ooov = transform_eri(ERI, Co, Co, Co, Cv)
    ooov = T.(ooov)

    oovv = transform_eri(ERI, Co, Co, Cv, Cv)
    oovv = T.(oovv)

    ovov = transform_eri(ERI, Co, Cv, Co, Cv)
    ovov = T.(ovov)

    ovvv = transform_eri(ERI, Co, Cv, Cv, Cv)
    ovvv = T.(ovvv)

    vvvv = transform_eri(ERI, Cv, Cv, Cv, Cv)
    vvvv = T.(vvvv)


    # Get density matrix

    D = Fermi.contract(C[:,1:ndocc], C[:,1:ndocc], "um", "vm")
    
    # Get AO Fock Matrix

    F = aoint.T + aoint.V
    Fermi.contract!(F,D,aoint.ERI,1.0,1.0,2.0,"mn","rs","mnrs")
    Fermi.contract!(F,D,aoint.ERI,1.0,1.0,-1.0,"mn","rs","mrns")

    # Get MO Fock Matrices

    oo = T.(transform_fock(F, Co, Co))
    ov = T.(transform_fock(F, Co, Cv))
    vv = T.(transform_fock(F, Cv, Cv))

    return PhysRestrictedMOIntegrals{T}(oooo, ooov, oovv, ovov, ovvv, vvvv, oo, ov, vv)
end

function transform_fock(F::Array{Float64,2}, C1::Array{Float64,2}, C2::Array{Float64,2})

    nmo,p = size(C1)
    _,q   = size(C2)

    Q1 = zeros(p,nmo)
    Fermi.contract!(Q1,C1,F,"pv","up","uv")

    Q2 = zeros(p,q)
    Fermi.contract!(Q2,C2,Q1,"pq","vq","pv")

    return Q2
end

function transform_eri(ERI::Fermi.MemTensor, C1::Array{Float64,2}, C2::Array{Float64,2}, C3::Array{Float64,2}, C4::Array{Float64,2})

    nmo,i = size(C1)
    _,a   = size(C2)
    _,j   = size(C3)
    _,b   = size(C4)

    Q1 = zeros(i,nmo,nmo,nmo)
    Fermi.contract!(Q1,C1,ERI,"ivls","ui","uvls")

    Q2 = zeros(i,a,nmo,nmo)
    Fermi.contract!(Q2,C2,Q1,"ials","va","ivls")
    Q1 = nothing

    Q3 = zeros(i,a,j,nmo)
    Fermi.contract!(Q3,C3,Q2,"iajs","lj","ials")
    Q2 = nothing

    Q4 = zeros(i,a,j,b)
    Fermi.contract!(Q4,C4,Q3,"iajb","sb","iajs")

    return Q4
end

end #module
