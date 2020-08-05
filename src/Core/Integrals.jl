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

import Base.getindex
import Base.setindex!

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
    bname::String
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

function aokinetic(molecule::Molecule, basis::String)#, interconnect::Fermi.Environments.No_IC,
                                                     #   communicator::Fermi.Environments.NoCommunicator,
                                                     #   accelerator::Fermi.Environments.NoAccelerator)

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
    T_engine = Lints.KineticEngine(nprim,l)
    sz = Lints.nao(bas)
    T = zeros(sz,sz)
    Lints.make_2D(T,T_engine,bas)
    T,bas
end
function aooverlap(molecule::Molecule, basis::String)#, interconnect::Fermi.Environments.No_IC,
                                                     #   communicator::Fermi.Environments.NoCommunicator,
                                                     #   accelerator::Fermi.Environments.NoAccelerator)

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
    sz = Lints.nao(bas)
    S = zeros(sz,sz)
    Lints.make_2D(S,S_engine,bas)
    S,bas
end
function aonuclear(molecule::Molecule, basis::String)#, interconnect::Fermi.Environments.No_IC,
                                                     #   communicator::Fermi.Environments.NoCommunicator,
                                                     #   accelerator::Fermi.Environments.NoAccelerator)

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
    V_engine = Lints.NuclearEngine(nprim,l,mol)
    sz = Lints.nao(bas)
    V = zeros(sz,sz)
    Lints.make_2D(V,V_engine,bas)
    V,bas
end
function aoeri(molecule::Molecule, basis::String)#, interconnect::Fermi.Environments.No_IC,
                                                 #       communicator::Fermi.Environments.NoCommunicator,
                                                 #       accelerator::Fermi.Environments.NoAccelerator)

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
    I_engines = []
    sz = Lints.nao(bas)
    for i in 1:Threads.nthreads()
        push!(I_engines,Lints.ERIEngine(nprim,l))
    end
    I = zeros(sz,sz,sz,sz)
    Lints.make_ERI(I,I_engines,bas)
    I,bas
end

struct DFAOIntegrals{T} <: AbstractAOIntegrals where T <: AbstractFloat
    bname::String
    dfbname::String
    basis::Lints.BasisSetAllocated
    dfbasis::Lints.BasisSetAllocated
    S::Array{T,2}
    T::Array{T,2}
    V::Array{T,2}
    ERI::Array{T,3}
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

function dfaoeri(molecule::Molecule, bname::String,dfbname::String)
    open("/tmp/molfile.xyz","w") do molfile
        natom = length(molecule.atoms)
        write(molfile,"$natom\n\n")
        write(molfile,Fermi.Geometry.get_xyz(molecule))
    end

    #@assert lowercase(bname) in ["cc-pvdz",
    #                              "cc-pvtz",
    #                              "cc-pvqz",
    #                              "cc-pv5z",
    #                              "aug-cc-pvdz",
    #                              "aug-cc-pvtz",
    #                              "aug-cc-pvqz",
    #                              "aug-cc-pv5z"] "Only Dunning basis sets are supported for density fitting!"
    #dfbname = bname*"-RI"
    Lints.libint2_init()
    #@output "   • Lints started\n\n"
    #@output "   Basis set: {}\n" bname
    #@output "   DF Basis set: {}\n" dfbname 
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
    Jh = Array(Hermitian(J)^(-1/2)) #sometimes Jh becomes complex slightly if J is not ~~exactly~~ hermitian 💔
    B = zeros(dfsz,sz,sz)
    Fermi.contract!(B,Pqp,Jh,"Qpq","Pqp","PQ")
    #@output "Exiting Lints.\n\n"
    B,dfbas
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

struct IntegralHelper{T}
    cache::Dict{String,Array} 
    bname::Dict{String,String}
    mol::Molecule
    C::Dict{String,Array{Float64,2}}
    basis::Dict{String,Lints.BasisSetAllocated} 
    type::DataType
end

function IntegralHelper()
    IntegralHelper{Float64}()
end

function IntegralHelper{T}() where T <: AbstractFloat
    cache = Dict{String,Array}() 
    bname = Dict{String,String}()
    type = T
    bname["primary"] = Fermi.CurrentOptions["basis"]
    aux = Fermi.CurrentOptions["auxbasis"]
    if aux == "auto"
        aux_lookup = Dict{String,String}(
                                         "cc-pvdz" => "cc-pvdz-jkfit",
                                         "cc-pvtz" => "cc-pvtz-jkfit",
                                         "cc-pvqz" => "cc-pvqz-jkfit",
                                         "cc-pv5z" => "cc-pv5z-jkfit"
                                        )
        bname["aux"] = try
            aux_lookup[Fermi.CurrentOptions["basis"]]
        catch KeyError #if we haven't got it programmed, use a large DF basis by default
            "aug-cc-pvqz-ri"
        end
    else
        bname["aux"] = aux
    end
    mol = Molecule()
    C = Dict{String,Array{Float64,2}}()
    basis = Dict{String,Lints.BasisSetAllocated}()
    IntegralHelper{T}(cache,bname,mol,C,basis,type)
end

function aux_ri!(I::IntegralHelper,ri=Fermi.CurrentOptions["rifit"])
    delete!(I.cache,"B") #clear out old aux basis
    if ri == "auto"
        aux_lookup = Dict{String,String}(
                                         "cc-pvdz" => "cc-pvdz-ri",
                                         "cc-pvtz" => "cc-pvtz-ri",
                                         "cc-pvqz" => "cc-pvqz-ri",
                                         "cc-pv5z" => "cc-pv5z-ri"
                                        )
        I.bname["aux"] = try
            aux_lookup[Fermi.CurrentOptions["basis"]]
        catch KeyError
            "aug-cc-pvqz-ri" # default to large DF basis
        end
    else
        I.bname["aux"] = ri
    end
end

"""
     O               -> occ,alpha,ERI
     o               -> occ,beta,ERI
     V               -> vir,alpha,ERI
     v               -> vir,beta,ERI
     P               -> all,alpha,ERI
     p               -> all,beta,ERI
     B               -> DF-ERI
     μ               -> ao,ERI
     Ω               -> NO,occ,alpha
     ω               -> NO,occ,beta
     U               -> NO,vir,alpha
     u               -> NO,vir,beta
     S               -> AO,overlap
     T               -> AO,kinetic
     V               -> AO,nuclear
    
     examples

"""
function getindex(I::IntegralHelper,entry::String)
    try
        I.cache[entry]
    catch KeyError
        compute!(I,entry)
        I.cache[entry]
    end
end

function setindex!(I::IntegralHelper,val,entry)
    I.cache[entry] = val
end

function compute!(I::IntegralHelper,entry::String)
    if entry == "μ" #AO basis eri. 4 index
        I.cache["μ"],I.basis["primary"] = aoeri(I.mol,I.bname["primary"])  

    elseif entry == "B" #AO basis eri. 3 index. DF
        I.cache["B"],I.basis["aux"] = dfaoeri(I.mol,I.bname["primary"],I.bname["aux"])

    elseif entry == "S" #AO basis overlap
        I.cache["S"],_ = aooverlap(I.mol,I.bname["primary"])

    elseif entry == "T" #AO basis kinetic
        I.cache["T"],_ = aokinetic(I.mol,I.bname["primary"])

    elseif entry == "V" #AO basis nuclear
        I.cache["V"],_ = aonuclear(I.mol,I.bname["primary"])

    elseif entry == "F"
        D = Fermi.contract(I.C["O"],I.C["O"],"um","vm")
        F = zeros(I.type,size(I["S"]))
        Fermi.HartreeFock.build_fock!(F,
                                                     I["T"] + I["V"],
                                                     D,
                                                     I["μ"])
        I.cache["F'"] = F
    elseif entry == "F'"
        D = Fermi.contract(I.C["O"],I.C["O"],"um","vm")
        F = zeros(I.type,size(I["S"]))
        Fermi.HartreeFock.build_fock!(F,
                                                     I["T"] + I["V"],
                                                     D,
                                                     I["B"])
        I.cache["F'"] = F
    elseif 'B' in entry #MO basis eri. 3 index. DF
        aoint = I["B"]
        C1 = try 
            I.C[entry[2:2]]
        catch
            error("")
        end
        C2 = try
            I.C[entry[3:3]]
        catch
            error("")
        end
        C = []
        push!(C,C1)
        push!(C,C2)
        drop_occ = Fermi.CurrentOptions["drop_occ"]
        drop_vir = Fermi.CurrentOptions["drop_vir"]

        for i in eachindex(C)
            C[i] == I.C["O"] ? C[i] = C[i][:, (1+drop_occ):end] : nothing
            C[i] == I.C["V"] ? C[i] = C[i][:, 1:end-drop_vir] : nothing
        end
        I.cache[entry] = transform_eri(aoint, C...)

    elseif 'F' in entry
        df = '\'' in entry
        if df
            F =  I["F'"]
        else
            F = I["F"]
        end
        C1 = try 
            I.C[entry[2:2]]
        catch
            error("")
        end
        C2 = try
            I.C[entry[3:3]]
        catch
            error("")
        end
        C = []
        push!(C,C1)
        push!(C,C2)
        drop_occ = Fermi.CurrentOptions["drop_occ"]
        drop_vir = Fermi.CurrentOptions["drop_vir"]

        for i in eachindex(C)
            C[i] == I.C["O"] ? C[i] = C[i][:, (1+drop_occ):end] : nothing
            C[i] == I.C["V"] ? C[i] = C[i][:, 1:end-drop_vir] : nothing
        end
        I[entry] = transform_fock(F,C...)

    else 
        aoint  = I["μ"] 
        if ('<' in entry || '(' in entry)
            offset = 1
        else
            offset = 0
        end
        if '(' in entry
            notation = "chem"
        else
            notation = "phys"
        end
        C = []
        C1 = try 
            I.C[entry[1+offset:1+offset]]
        catch
            error("")
        end
        C2 = try
            I.C[entry[2+offset:2+offset]]
        catch
            error("")
        end
        C3 = try
            I.C[entry[3+offset:3+offset]]
        catch
            error("")
        end
        C4 = try
            I.C[entry[4+offset:4+offset]]
        catch
            error("")
        end
        push!(C,C1)
        push!(C,C3)
        push!(C,C2)
        push!(C,C4)
        drop_occ = Fermi.CurrentOptions["drop_occ"]
        drop_vir = Fermi.CurrentOptions["drop_vir"]

        for i in eachindex(C)
            C[i] == I.C["O"] ? C[i] = C[i][:, (1+drop_occ):end] : nothing
            C[i] == I.C["V"] ? C[i] = C[i][:, 1:end-drop_vir] : nothing
        end
        temp = transform_eri(aoint, C...)
        if notation == "phys"
            I.cache[entry] = permutedims(temp,(1,3,2,4))
        else
            I.cache[entry] = temp
        end
    end
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

function transform_eri(ERI::Array{T,4}, C1::Array{Float64,2}, C2::Array{Float64,2}, C3::Array{Float64,2}, C4::Array{Float64,2}) where T <: AbstractFloat

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
    Q4
end

function transform_eri(ERI::Array{T,3},C1::Array{T,2},C2::Array{T,2}) where T <: AbstractFloat
    naux = size(ERI,1)
    nao,i = size(C1)
    _,a = size(C2)
    Bin = zeros(naux,i,nao)
    Fermi.contract!(Bin,C1,ERI,"Qiv","ui","Quv")
    Bia = zeros(naux,i,a)
    Fermi.contract!(Bia,C2,Bin,"Qia","va","Qiv")
    Bia
end

end #module
