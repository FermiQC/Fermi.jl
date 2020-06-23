"""
Module for storing and handling reference wavefunctions.

## structs

    Wfn -> holds info for disk based and in core computations.
    DirectWfn -> holds info for integral direct computations.

"""
module Wavefunction

using Fermi.DiskTensors
using Fermi.Transformation
using Fermi


export Wfn
export DirectWfn

"""
	Wfn
Data structure for storing integrals, MO coefficients, and misc information about
a reference (HF) wavefunction.

## Fields
energy::T energy of the reference wavefunction

vnuc::T nuclear repulsion

nalpha::Int number of alpha electrons

nbeta::Int number of beta electrons

nvira::Int number of virtual functions of alpha spin

nvirb::Int number of virtual functions of beta spin

nmo::Int number of molecular orbitals

unrestricted::Bool whether or not the alpha and beta spatial extents are required to
be the same.

Ca::Array{T,2} AO->MO coefficients for alpha MO's

Cb::Array{T,2} AO->MO coefficients for beta MO's

hao::Array{T,2} AO basis core hamiltonian (kinetic + potential)

epsa::Array{T,1} orbital eigenvalues for alpha MO's

epsb::Array{T,1} orbtial eigenvalues for beta MO's

ao_eri::Union{Array{T,4},DiskFourTensor} AO basis TEI

pqrs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case AAAA)

pQrS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case ABAB)

pQRs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case ABBA)

PQRS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BBBB)

PqRs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BABA)

PqrS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BAAB)
"""
struct Wfn{T}
    energy::T
    vnuc::T
    nalpha::Int
    nbeta::Int
    nvira::Int
    nvirb::Int
    nmo::Int
    unrestricted::Bool
    basis
    molecule
    Ca::Array{T,2} #AO->MO coefficients
    Cb::Array{T,2} #AO->MO coefficients
    Cao::Array{T,2}
    Cav::Array{T,2}
    Cbo::Array{T,2}
    Cbv::Array{T,2}
    hao::Array{T,2} #Core hamiltonian
    epsa::Array{T,1} #orbital eigenvalues
    epsb::Array{T,1} #orbital eigenvalues
    ao_eri::Union{Array{T,4},DiskFourTensor} #AO basis electron repulsion integrals
end

function Wfn(rhfwfn::Fermi.HartreeFock.RHF.RHFWfn)
    Wfn{Float64}(rhfwfn)
end

function Wfn{T}(wfn::Fermi.HartreeFock.RHF.RHFWfn; unrestricted::Bool=false, diskbased::Bool=false, name::String = "default", df::Bool=false) where T
    dt = T
    energy = wfn.energy[1]
    dummy2 = Array{dt}(undef, 0, 0) #placeholder 2D array
    dummy4 = Array{dt}(undef, 0, 0, 0, 0) #placeholder 4D array
    vnuc = wfn.vnuc
    basis = wfn.basis
    molecule = wfn.molecule
    nbf = size(wfn.C,1)
    nocca = fld(wfn.nelec,2)
    nvira = nbf - nocca
    noccb = fld(wfn.nelec,2)
    nvirb = nbf - noccb
    epsa = convert(Array{dt,1}, wfn.eps)
    epsb = convert(Array{dt,1}, wfn.eps)
    Cao = convert(Array{dt,2}, wfn.C[:,1:nocca])
    Cav = convert(Array{dt,2}, wfn.C[:,nocca+1:end])
    Cbo = convert(Array{dt,2}, wfn.C[:,1:nocca])
    Cbv = convert(Array{dt,2}, wfn.C[:,nocca+1:end])
    Ca = convert(Array{dt,2}, wfn.C)
    Cb = convert(Array{dt,2}, wfn.C)
    hao = convert(Array{dt,2}, wfn.H) #core hamiltonian in AO
     
    # create the Wfn object and return it!
    owfn = Wfn{dt}(
        energy,
        vnuc,
        nocca,
        noccb,
        nvira,
        nvirb,
        nbf,
        unrestricted,
        basis,
        molecule,
        Ca,
        Cb,
        Cao,
        Cav,
        Cbo,
        Cbv,
        hao,
        epsa,
        epsb,
        wfn.I
    )
    return owfn
end

end # Module
