"""
    Fermi.Integrals

Module to compute integrals using Lints.jl
"""
module Integrals

using Fermi
using Fermi.Output
using Fermi.Geometry: Molecule
using LinearAlgebra
using TensorOperations

import Base: getindex, setindex!, delete!

include("Lints.jl")

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
    mol::Molecule
    basis::String
    aux::String
    cache::Dict{String,FermiMDArray{T}} 
    normalize::Bool
end

function IntegralHelper(x...)
    IntegralHelper{Float64}(x...)
end

function IntegralHelper{T}() where T <: AbstractFloat
    mol = Molecule()
    basis = Fermi.CurrentOptions["basis"]
    aux = Fermi.CurrentOptions["jkfit"]
    IntegralHelper{T}(mol, basis, aux)
end

function IntegralHelper{T}(mol::Molecule) where T <: AbstractFloat
    basis = Fermi.CurrentOptions["basis"]
    aux = Fermi.CurrentOptions["jkfit"]
    IntegralHelper{T}(mol, basis, aux)
end

function IntegralHelper{T}(mol::Molecule, basis::String) where T <: AbstractFloat
    aux = Fermi.CurrentOptions["jkfit"]
    IntegralHelper{T}(mol, basis, aux)
end

function IntegralHelper{T}(mol::Molecule, basis::String, aux::String) where T <: AbstractFloat
    if aux == "auto"
        aux_lookup = Dict{String,String}(
                                         "cc-pvdz" => "cc-pvdz-jkfit",
                                         "cc-pvtz" => "cc-pvtz-jkfit",
                                         "cc-pvqz" => "cc-pvqz-jkfit",
                                         "cc-pv5z" => "cc-pv5z-jkfit"
                                        )
        aux = haskey(aux_lookup, basis) ? aux_lookup[basis] : "aug-cc-pvqz-rifit"
    end
    cache = Dict{String, FermiMDArray{T}}() 
    lbasis = Dict{String,Lints.BasisSetAllocated}()
    IntegralHelper{T}(mol, basis, aux, cache, false)
end

# Clears cache and change normalize key
function normalize!(I::IntegralHelper,normalize::Bool)
    if I.normalize != normalize
        I.normalize = normalize
        for entry in keys(I.cache)
            delete!(I.cache, entry)
        end
    end
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
            "aug-cc-pvqz-rifit" # default to large DF basis
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
function aux_jk!(I::IntegralHelper, jk = Fermi.CurrentOptions["jkfit"])
    delete!(I,"B") #clear out old aux basis
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
            "aug-cc-pvqz-rifit" # default to large DF basis
        end
    else
        I.aux = jk
    end
end

"""
    getindex(I::IntegralHelper,entry::String)

Called when `helper["foo"]` syntax is used. If the requested entry already
exists, simply return the entry. If not, compute the requested entry.
"""
function getindex(I::IntegralHelper,entry::String)
    if haskey(I.cache, entry)
        return I.cache[entry]
    else
        compute!(I, entry)
        return I.cache[entry]
    end
end

function delete!(I::IntegralHelper, key)
    delete!(I.cache, key)
    GC.gc()
end

function compute!(I::IntegralHelper, entry::String)

    if entry == "S" #AO basis overlap
        I.cache["S"] = FermiMDArray(ao_overlap(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "T" #AO basis kinetic
        I.cache["T"] = FermiMDArray(ao_kinetic(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "V" #AO basis nuclear
        I.cache["V"] = FermiMDArray(ao_nuclear(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "ERI" 
        I.cache["ERI"] = FermiMDArray(ao_eri(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "DFERI"
        I.cache["DFERI"] = FermiMDArray(df_ao_eri(I.mol, I.basis, I.aux, normalize = I.normalize))

    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

#function transform_fock(F::Array{Float64,2}, O1::O, O2::O) where O <: AbstractOrbitals
#    C1 = hcat([orb.C for orb in O1.orbs]...)
#    C2 = hcat([orb.C for orb in O2.orbs]...)
#    transform_fock(F,C1,C2)
#end

#function transform_fock(F::Array{Float64,2}, C1::Array{Float64,2}, C2::Array{Float64,2})
#
#    nmo,p = size(C1)
#    _,q   = size(C2)
#
#    Q1 = zeros(p,nmo)
#    Fermi.contract!(Q1,C1,F,"pv","up","uv")
#
#    Q2 = zeros(p,q)
#    Fermi.contract!(Q2,C2,Q1,"pq","vq","pv")
#
#    return Q2
#end
#
#function transform_eri(ERI::Array{T,4}, O1::O, O2::O, O3::O, O4::O) where { O <: AbstractOrbitals,
#                                                                            T <: AbstractFloat }
#    C1 = hcat([orb.C for orb in O1.orbs]...)
#    C2 = hcat([orb.C for orb in O2.orbs]...)
#    C3 = hcat([orb.C for orb in O3.orbs]...)
#    C4 = hcat([orb.C for orb in O4.orbs]...)
#    transform_eri(ERI,C1,C2,C3,C4)
#end
#
#function transform_eri(ERI::Array{T,3}, O1::O, O2::O) where { O <: AbstractOrbitals,
#                                                              T <: AbstractFloat }
#    C1 = hcat([orb.C for orb in O1.orbs]...)
#    C2 = hcat([orb.C for orb in O2.orbs]...)
#    transform_eri(ERI,C1,C2)
#end
#
#function transform_eri(ERI::Array{T,4}, C1::Array{Float64,2}, C2::Array{Float64,2}, C3::Array{Float64,2}, C4::Array{Float64,2}) where T <: AbstractFloat
#
#    nmo,i = size(C1)
#    _,a   = size(C2)
#    _,j   = size(C3)
#    _,b   = size(C4)
#
#    C1 = T.(C1)
#    C2 = T.(C2)
#    C3 = T.(C3)
#    C4 = T.(C4)
#    Q1 = zeros(T,i,nmo,nmo,nmo)
#    Fermi.contract!(Q1,C1,ERI,"ivls","ui","uvls")
#
#    Q2 = zeros(T,i,a,nmo,nmo)
#    Fermi.contract!(Q2,C2,Q1,"ials","va","ivls")
#    Q1 = nothing
#
#    Q3 = zeros(T,i,a,j,nmo)
#    Fermi.contract!(Q3,C3,Q2,"iajs","lj","ials")
#    Q2 = nothing
#
#    Q4 = zeros(T,i,a,j,b)
#    Fermi.contract!(Q4,C4,Q3,"iajb","sb","iajs")
#    Q4
#end
#
#function transform_eri(ERI::Array{T,3},C1::Array{Float64,2},C2::Array{Float64,2}) where T <: AbstractFloat
#    C1 = T.(C1)
#    C2 = T.(C2)
#    naux = size(ERI,1)
#    nao,i = size(C1)
#    _,a = size(C2)
#    Bin = zeros(T,naux,i,nao)
#    Fermi.contract!(Bin,C1,ERI,"Qiv","ui","Quv")
#    Bia = zeros(T,naux,i,a)
#    Fermi.contract!(Bia,C2,Bin,"Qia","va","Qiv")
#    Bia
#end

end #module
