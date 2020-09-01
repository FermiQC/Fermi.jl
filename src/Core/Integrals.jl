"""
    Fermi.Integrals

Module to compute integrals using Lints.jl
"""
module Integrals

using Fermi
using Fermi.Output
using Fermi.Geometry: Molecule
using Fermi.Orbitals: AbstractOrbitals, OrbDict
using Lints
using LinearAlgebra
using TensorOperations

import Base.getindex
import Base.setindex!

"""
    IntegralHelper{T}

Structure to assist with computing and storing integrals. 
Accesss like a dictionary e.g.,
    ints["B"]
to get DF ERI integrals. A number of characters can be added to denote orbitals in various
bases, such as MO or NO.
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

# Fields
    cache::Dict{String,Array}                        previously computed integrals
    bname::Dict{String,String}                       basis set names for various purposes
    mol::Molecule                                    attached Molecule object.
    orbs::Dict{String,O} where O < AbstractOrbitals  orbitals of various kinds
    basis::Dict{String,Lints.BasisSetAllocated}      basis sets for various purposes
"""
mutable struct IntegralHelper{T}
    cache::Dict{String,Array} 
    bname::Dict{String,String}
    mol::Molecule
    orbs::Fermi.Orbitals.OrbDict
    basis::Dict{String,Lints.BasisSetAllocated} 
    type::DataType
end

function IntegralHelper(;mol=Molecule(),orbs=Fermi.Orbitals.OrbDict())
    IntegralHelper{Float64}(;mol=mol,orbs=orbs)
end

function IntegralHelper{T}(;mol=Molecule(),orbs=Fermi.Orbitals.OrbDict()) where T <: AbstractFloat
    cache = Dict{String,Array}() 
    bname = Dict{String,String}()
    type = T
    bname["primary"] = Fermi.CurrentOptions["basis"]
    aux = Fermi.CurrentOptions["jkfit"]
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
            "aug-cc-pvqz-rifit"
        end
    else
        bname["aux"] = aux
    end
    #mol = Molecule()
    #orbs = Fermi.Orbitals.OrbDict()
    basis = Dict{String,Lints.BasisSetAllocated}()
    IntegralHelper{T}(cache,bname,mol,orbs,basis,type)
end


"""
    aokinetic(molecule::Fermi.Molecule, basis::String)

Computes AO basis kinetic energy ⟨μ|T̂|ν⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["T"]

where `helper` is bound to the desired molecule and basis set.
"""
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
    Lints.libint2_finalize()
    T_engine = nothing
    GC.gc()
    T,bas
end
"""
    aooverlap(molecule::Fermi.Molecule, basis::String)

Computes AO basis overlap ⟨p|q⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["S"]

where `helper` is bound to the desired molecule and basis set.
"""
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
    Lints.libint2_finalize()
    S_engine = nothing
    GC.gc()
    S,bas
end
function aodipole(molecule::Molecule, basis::String)#, interconnect::Fermi.Environments.No_IC,
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
    S_engine = Lints.DipoleEngine(nprim,l)
    sz = Lints.nao(bas)
    S = zeros(sz,sz)
    Lints.make_2D(S,S_engine,bas)
    Lints.libint2_finalize()
    S_engine = nothing
    GC.gc()

    S,bas
end
"""
    aooverlap(molecule::Fermi.Molecule, basis::String)

Computes AO basis overlap ⟨μ|V̂|ν⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["V"]

where `helper` is bound to the desired molecule and basis set.
"""
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
    Lints.libint2_finalize()
    V_engine = nothing
    GC.gc()
    V,bas
end
"""
    aoeri(molecule::Fermi.Molecule, basis::String)

Computes AO basis electron repulsion integrals ⟨μν|Ô₂|ρσ⟩ integrals for the given basis and molecule.
Can be accessed at a higher level by calling
    
    helper["μ"]

where `helper` is bound to the desired molecule and basis set.
"""
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
    Lints.libint2_finalize()
    I_engines = nothing
    GC.gc()
    I,bas
end

"""
    dfaoeri(molecule::Fermi.Molecule, basis::String)

Computes AO basis density fitted electron repulsion integrals ⟨μν|Ô₂|P⟩J(p,q)^-1/2 integrals for the given basis and molecule.
Note that the returned integrals DO NOT need to be combined with the Coulomb metric J(P,Q). In common notation, this is B(Q,μ,ν).
Can be accessed at a higher level by calling
    
    helper["B"]

where `helper` is bound to the desired molecule and basis set.
"""
function dfaoeri(molecule::Molecule, bname::String,dfbname::String)
    open("/tmp/molfile.xyz","w") do molfile
        natom = length(molecule.atoms)
        write(molfile,"$natom\n\n")
        write(molfile,Fermi.Geometry.get_xyz(molecule))
    end

    Lints.libint2_init()
    mol = Lints.Molecule("/tmp/molfile.xyz")
    bas = Lints.BasisSet(bname, mol)
    dfbas = Lints.BasisSet(dfbname,mol)

    nprim = max(Lints.max_nprim(bas),Lints.max_nprim(dfbas))
    l = max(Lints.max_l(bas),Lints.max_l(dfbas))

    eri_engines = Any[Lints.DFEngine(nprim,l) for i=1:Threads.nthreads()]
    sz = Lints.nao(bas)
    dfsz = Lints.nao(dfbas)
    J = zeros(dfsz,dfsz)
    Pqp = zeros(dfsz,sz,sz)
    Lints.make_b(Pqp,eri_engines,bas,dfbas)
    Lints.make_j(J,eri_engines[1],dfbas)
    Jh = Array(Hermitian(J)^(-1/2)) #sometimes Jh becomes complex slightly if J is not ~~exactly~~ hermitian 💔
    B = zeros(dfsz,sz,sz)
    Lints.libint2_finalize()
    for i=1:Threads.nthreads()
        eri_engines[i] = nothing
    end
    eri_engines = nothing
    GC.gc()
    Fermi.contract!(B,Pqp,Jh,"Qpq","Pqp","PQ")
    B,dfbas
end

"""
    aux_ri!(I::IntegralHelper, ri=Fermi.CurrentOptions["rifit"])

Clears the integral cache and switches auxiliary DF integrals to use the
current RI fitting basis set. Used between DF-RHF and DF-post HF.
"""
function aux_ri!(I::IntegralHelper,ri=Fermi.CurrentOptions["rifit"])
    delete!(I.cache,"B") #clear out old aux basis
    GC.gc()
    if ri == "auto"
        aux_lookup = Dict{String,String}(
                                         "cc-pvdz" => "cc-pvdz-rifit",
                                         "cc-pvtz" => "cc-pvtz-rifit",
                                         "cc-pvqz" => "cc-pvqz-rifit",
                                         "cc-pv5z" => "cc-pv5z-rifit"
                                        )
        I.bname["aux"] = try
            aux_lookup[Fermi.CurrentOptions["basis"]]
        catch KeyError
            "augmentation-cc-pvqz-rifit" # default to large DF basis
        end
    else
        I.bname["aux"] = ri
    end
end
"""
    aux_ri!(I::IntegralHelper, jk=Fermi.CurrentOptions["jkfit"])

Clears the integral cache and switches auxiliary DF integrals to use the
current JK fitting basis set. Used to ensure JK integrals are used in
DF-RHF.
"""
function aux_jk!(I::IntegralHelper,jk=Fermi.CurrentOptions["jkfit"])
    delete!(I.cache,"B") #clear out old aux basis
    GC.gc()
    if jk == "auto"
        aux_lookup = Dict{String,String}(
                                         "cc-pvdz" => "cc-pvdz-jkfit",
                                         "cc-pvtz" => "cc-pvtz-jkfit",
                                         "cc-pvqz" => "cc-pvqz-jkfit",
                                         "cc-pv5z" => "cc-pv5z-jkfit"
                                        )
        I.bname["aux"] = try
            aux_lookup[Fermi.CurrentOptions["basis"]]
        catch KeyError
            "augmentation-cc-pvqz-rifit" # default to large DF basis
        end
    else
        I.bname["aux"] = jk
    end
end

"""
    getindex(I::IntegralHelper,entry::String)

Called when `helper["foo"]` syntax is used. If the requested entry already
exists, simply return the entry. If not, compute the requested entry.
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
    o = I.orbs
    if entry == "μ" #AO basis eri. 4 index
        I.cache["μ"],_ = aoeri(I.mol,I.bname["primary"])  

    elseif entry == "B" #AO basis eri. 3 index. DF
        I.cache["B"],_ = dfaoeri(I.mol,I.bname["primary"],I.bname["aux"])

    elseif entry == "S" #AO basis overlap
        I.cache["S"],_ = aooverlap(I.mol,I.bname["primary"])

    elseif entry == "T" #AO basis kinetic
        I.cache["T"],_ = aokinetic(I.mol,I.bname["primary"])

    elseif entry == "V" #AO basis nuclear
        I.cache["V"],_ = aonuclear(I.mol,I.bname["primary"])

    elseif entry == "F"
        D = Fermi.contract(I.C["[O]"],I.C["[O]"],"um","vm")
        F = zeros(I.type,size(I["S"]))
        Fermi.HartreeFock.build_fock!(F,
                                                     I["T"] + I["V"],
                                                     D,
                                                     I["μ"])
        I.cache["F'"] = F
    elseif entry == "F'"
        D = Fermi.contract(I.C["[O]"],I.C["[O]"],"um","vm")
        F = zeros(I.type,size(I["S"]))
        Fermi.HartreeFock.build_fock!(F,
                                                     I["T"] + I["V"],
                                                     D,
                                                     I["B"])
        I.cache["F'"] = F
    elseif 'B' in entry #MO basis eri. 3 index. DF
        aoint = I["B"]
        C1 = try 
            o[entry[2:2]]
        catch KeyError
            error("")
        end
        C2 = try
            o[entry[3:3]]
        catch KeyError
            error("")
        end
        C = []
        push!(C,C1)
        push!(C,C2)
        I.cache[entry] = transform_eri(aoint, C...)

    elseif 'F' in entry
        df = '\'' in entry
        if df
            F =  I["F'"]
        else
            F = I["F"]
        end
        C1 = try 
            o[entry[2:2]]
        catch
            error("")
        end
        C2 = try
            o[entry[3:3]]
        catch
            error("")
        end
        C = []
        push!(C,C1)
        push!(C,C2)
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
            o[entry[1+offset:1+offset]]
        catch
            error("")
        end
        C2 = try
            o[entry[2+offset:2+offset]]
        catch
            error("")
        end
        C3 = try
            o[entry[3+offset:3+offset]]
        catch
            error("")
        end
        C4 = try
            o[entry[4+offset:4+offset]]
        catch
            error("")
        end
        push!(C,C1)
        push!(C,C3)
        push!(C,C2)
        push!(C,C4)
        temp = transform_eri(aoint, C...)
        if notation == "phys"
            I.cache[entry] = permutedims(temp,(1,3,2,4))
        else
            I.cache[entry] = temp
        end
    end
end

function transform_fock(F::Array{Float64,2}, O1::O, O2::O) where O <: AbstractOrbitals
    C1 = hcat([orb.C for orb in O1.orbs]...)
    C2 = hcat([orb.C for orb in O2.orbs]...)
    transform_fock(F,C1,C2)
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

function transform_eri(ERI::Array{T,4}, O1::O, O2::O, O3::O, O4::O) where { O <: AbstractOrbitals,
                                                                            T <: AbstractFloat }
    C1 = hcat([orb.C for orb in O1.orbs]...)
    C2 = hcat([orb.C for orb in O2.orbs]...)
    C3 = hcat([orb.C for orb in O3.orbs]...)
    C4 = hcat([orb.C for orb in O4.orbs]...)
    transform_eri(ERI,C1,C2,C3,C4)
end

function transform_eri(ERI::Array{T,3}, O1::O, O2::O) where { O <: AbstractOrbitals,
                                                              T <: AbstractFloat }
    C1 = hcat([orb.C for orb in O1.orbs]...)
    C2 = hcat([orb.C for orb in O2.orbs]...)
    transform_eri(ERI,C1,C2)
end

function transform_eri(ERI::Array{T,4}, C1::Array{Float64,2}, C2::Array{Float64,2}, C3::Array{Float64,2}, C4::Array{Float64,2}) where T <: AbstractFloat

    nmo,i = size(C1)
    _,a   = size(C2)
    _,j   = size(C3)
    _,b   = size(C4)

    C1 = T.(C1)
    C2 = T.(C2)
    C3 = T.(C3)
    C4 = T.(C4)
    Q1 = zeros(T,i,nmo,nmo,nmo)
    Fermi.contract!(Q1,C1,ERI,"ivls","ui","uvls")

    Q2 = zeros(T,i,a,nmo,nmo)
    Fermi.contract!(Q2,C2,Q1,"ials","va","ivls")
    Q1 = nothing

    Q3 = zeros(T,i,a,j,nmo)
    Fermi.contract!(Q3,C3,Q2,"iajs","lj","ials")
    Q2 = nothing

    Q4 = zeros(T,i,a,j,b)
    Fermi.contract!(Q4,C4,Q3,"iajb","sb","iajs")
    Q4
end

function transform_eri(ERI::Array{T,3},C1::Array{Float64,2},C2::Array{Float64,2}) where T <: AbstractFloat
    C1 = T.(C1)
    C2 = T.(C2)
    naux = size(ERI,1)
    nao,i = size(C1)
    _,a = size(C2)
    Bin = zeros(T,naux,i,nao)
    Fermi.contract!(Bin,C1,ERI,"Qiv","ui","Quv")
    Bia = zeros(T,naux,i,a)
    Fermi.contract!(Bia,C2,Bin,"Qia","va","Qiv")
    Bia
end

end #module
