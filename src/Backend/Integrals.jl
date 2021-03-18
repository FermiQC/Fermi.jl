"""
    Fermi.Integrals

Module to compute integrals using Lints.jl
"""
module Integrals

using Fermi
using Fermi.Geometry: Molecule
using LinearAlgebra
using TensorOperations

import Base: getindex, setindex!, delete!

export IntegralHelper
export delete!
export ao_to_mo_eri
export ao_to_mo_rieri
export ao_to_mo_eri!
export ao_to_mo_rieri!

include("Lints.jl")

"""
    IntegralHelper{T}

Structure to assist with computing and storing integrals. 
Accesss like a dictionary e.g.,
    ints["S"]

A key is associated with each type of integral

    "S"           -> AO overlap integral
    "T"           -> AO electron kinetic energy integral
    "V"           -> AO electron-nuclei attraction integral
    "ERI"         -> AO electron repulsion integral
    "JKERI"       -> AO JK density fitted electron repulsion integral
    "RIERI"       -> AO RI density fitted electron repulsion integral

# Fields
    mol                         Associated Fermi Molecule object
    basis                       Basis set used within the helper
    aux                         Auxiliar basis set used in density fitting
    cache                       Holds integrals already computed 
    notation                    `chem` or `phys` for ERI notation
    normalize                   Do normalize integrals? `true` or `false`
"""
mutable struct IntegralHelper{T}
    mol::Molecule
    basis::String
    auxjk::String
    auxri::String
    cache::Dict{String,FermiMDArray{T}} 
    notation::String
    normalize::Bool
end

function IntegralHelper(x...)
    IntegralHelper{Float64}(x...)
end

function IntegralHelper{T}() where T <: AbstractFloat
    mol = Molecule()
    basis = Fermi.Options.get("basis")
    auxjk = Fermi.Options.get("jkfit")
    auxri = Fermi.Options.get("rifit")
    IntegralHelper{T}(mol, basis, auxjk, auxri)
end

function IntegralHelper{T}(mol::Molecule) where T <: AbstractFloat
    basis = Fermi.Options.get("basis")
    auxjk = Fermi.Options.get("jkfit")
    auxri = Fermi.Options.get("rifit")
    IntegralHelper{T}(mol, basis, auxjk, auxri)
end

function IntegralHelper{T}(mol::Molecule, basis::String) where T <: AbstractFloat
    auxjk = Fermi.Options.get("jkfit")
    auxri = Fermi.Options.get("rifit")
    IntegralHelper{T}(mol, basis, auxjk, auxri)
end

function IntegralHelper{T}(mol::Molecule, basis::String, auxjk::String, auxri::String) where T <: AbstractFloat
    if auxjk == "auto"
        std_name = Regex("cc-pv.z")
        auxjk = occursin(std_name, basis) ? basis*"-jkfit" : "cc-pvqz-jkfit"
    end

    if auxri == "auto"
        std_name = Regex("cc-pv.z")
        auxri = occursin(std_name, basis) ? basis*"-rifit" : "cc-pvqz-rifit"
    end
    cache = Dict{String, FermiMDArray{T}}() 
    IntegralHelper{T}(mol, basis, auxjk, auxri, cache, "chem", false)
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

function delete!(I::IntegralHelper, keys...)
    for k in keys
        delete!(I.cache, k)
    end
    GC.gc()
end

function compute!(I::IntegralHelper{Float64}, entry::String)

    if entry == "S" #AO basis overlap
        I.cache["S"] = FermiMDArray(ao_overlap(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "T" #AO basis kinetic
        I.cache["T"] = FermiMDArray(ao_kinetic(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "V" #AO basis nuclear
        I.cache["V"] = FermiMDArray(ao_nuclear(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "ERI" 
        I.cache["ERI"] = FermiMDArray(ao_eri(I.mol, I.basis, normalize = I.normalize))

    elseif entry == "JKERI"
        I.cache["JKERI"] = FermiMDArray(df_ao_eri(I.mol, I.basis, I.auxjk, normalize = I.normalize))

    elseif entry == "RIERI"
        I.cache["RIERI"] = FermiMDArray(df_ao_eri(I.mol, I.basis, I.auxri, normalize = I.normalize))

    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute!(I::IntegralHelper{Float32}, entry::String)

    if entry == "S" #AO basis overlap
        I.cache["S"] = FermiMDArray(Float32.(ao_overlap(I.mol, I.basis, normalize = I.normalize)))

    elseif entry == "T" #AO basis kinetic
        I.cache["T"] = FermiMDArray(Float32.(ao_kinetic(I.mol, I.basis, normalize = I.normalize)))

    elseif entry == "V" #AO basis nuclear
        I.cache["V"] = FermiMDArray(Float32.(ao_nuclear(I.mol, I.basis, normalize = I.normalize)))

    elseif entry == "ERI" 
        I.cache["ERI"] = FermiMDArray(Float32.(ao_eri(I.mol, I.basis, normalize = I.normalize)))

    elseif entry == "JKERI"
        I.cache["JKERI"] = FermiMDArray(Float32.(df_ao_eri(I.mol, I.basis, I.auxjk, normalize = I.normalize)))

    elseif entry == "RIERI"
        I.cache["RIERI"] = FermiMDArray(Float32.(df_ao_eri(I.mol, I.basis, I.auxri, normalize = I.normalize)))

    else
        throw(Fermi.InvalidFermiOption("Invalid key for IntegralHelper: $(entry)."))
    end
end

function ao_to_mo_eri(I::IntegralHelper, C::AbstractArray{T,2}; name::String="MOERI") where T <: AbstractFloat

    if !(haskey(I.cache, "ERI"))
        output("No previous ERI found.")
        output("Computing ERI...")
        t = @elapsed I["ERI"]
        output("Done in {:5.5f}", t)
    end

    AOERI = I["ERI"]
    @tensoropt MOERI[p,q,r,s] :=  AOERI[μ, ν, ρ, σ]*C[μ, p]*C[ν, q]*C[ρ, r]*C[σ, s]

    I.cache[name] = MOERI
end

function ao_to_mo_eri(I::IntegralHelper, C1::AbstractArray{T,2}, C2::AbstractArray{T,2}, 
                     C3::AbstractArray{T,2}, C4::AbstractArray{T,2}; name::String="MOERI") where T <: AbstractFloat

    if !(haskey(I.cache, "ERI"))
        output("No previous ERI found.")
        output("Computing ERI...")
        t = @elapsed I["ERI"]
        output("Done in {:5.5f}", t)
    end

    AOERI = I["ERI"]
    @tensoropt MOERI[p,q,r,s] :=  AOERI[μ, ν, ρ, σ]*C1[μ, p]*C2[ν, q]*C3[ρ, r]*C4[σ, s]

    I.cache[name] = MOERI
end

function ao_to_mo_eri!(I::IntegralHelper, C...; name="MOERI") where T <: AbstractFloat

    ao_to_mo_eri(I, C..., name=name)
    delete!(I, "ERI")

    return I[name]
end

function ao_to_mo_rieri(I::IntegralHelper, C1::AbstractArray{T,2}, C2::AbstractArray{T,2}; name::String="MORIERI") where T <: AbstractFloat

    if !(haskey(I.cache, "RIERI"))
        output("No previous RI-ERI found.")
        output("Computing RI-ERI...")
        t = @elapsed I["RIERI"]
        output("Done in {:5.5f}", t)
    end

    AOERI = I["RIERI"]
    @tensoropt MOERI[P,p,q] :=  AOERI[P,μ, ν]*C1[μ, p]*C2[ν, q]

    I.cache[name] = MOERI
end

function ao_to_mo_rieri!(I::IntegralHelper, C1::AbstractArray{T,2}, C2::AbstractArray{T,2}; name::String="MORIERI") where T <: AbstractFloat

    ao_to_mo_rieri(I, C1, C2, name=name)
    delete!(I, "RIERI")

    return I[name]
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

#"""
#    aux_ri!(I::IntegralHelper, ri=Fermi.Options.get["rifit"])
#
#Clears the integral cache and switches auxiliary DF integrals to use the
#current RI fitting basis set. Used between DF-RHF and DF-post HF.
#"""
#function aux_ri!(I::IntegralHelper,ri=Fermi.Options.get("rifit"))
#    delete!(I.cache,"B") #clear out old aux basis
#    GC.gc()
#    if ri == "auto"
#        aux_lookup = Dict{String,String}(
#                                         "cc-pvdz" => "cc-pvdz-rifit",
#                                         "cc-pvtz" => "cc-pvtz-rifit",
#                                         "cc-pvqz" => "cc-pvqz-rifit",
#                                         "cc-pv5z" => "cc-pv5z-rifit"
#                                        )
#        I.bname["aux"] = try
#            aux_lookup[Fermi.Options.get("basis")]
#        catch KeyError
#            "aug-cc-pvqz-rifit" # default to large DF basis
#        end
#    else
#        I.bname["aux"] = ri
#    end
#end
#"""
#    aux_ri!(I::IntegralHelper, jk=Fermi.Options.get["jkfit"])
#
#Clears the integral cache and switches auxiliary DF integrals to use the
#current JK fitting basis set. Used to ensure JK integrals are used in
#DF-RHF.
#"""
#function aux_jk!(I::IntegralHelper, jk = Fermi.Options.get("jkfit"))
#    delete!(I,"B") #clear out old aux basis
#    if jk == "auto"
#        aux_lookup = Dict{String,String}(
#                                         "cc-pvdz" => "cc-pvdz-jkfit",
#                                         "cc-pvtz" => "cc-pvtz-jkfit",
#                                         "cc-pvqz" => "cc-pvqz-jkfit",
#                                         "cc-pv5z" => "cc-pv5z-jkfit"
#                                        )
#        I.bname["aux"] = try
#            aux_lookup[Fermi.Options.get("basis")]
#        catch KeyError
#            "aug-cc-pvqz-rifit" # default to large DF basis
#        end
#    else
#        I.aux = jk
#    end
#end

end #module
