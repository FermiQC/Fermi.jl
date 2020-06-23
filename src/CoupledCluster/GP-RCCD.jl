"""
    Fermi.CoupledCluster.GPRCCD

module for performing GPU accelerated RCCD computations

## Functions

    Fermi.CoupledCluster.GPRCCD.do_rccd
"""
module GPRCCD
using Fermi.Wavefunction
#using CuArrays
using TensorOperations
using Fermi.Transformation
include("Denominators.jl")
#println(CuArrays.has_cutensor())
export do_rccd
"""
    do_rccd

Applies RCCD equations to the input Wfn object. Disk based versus in-core 
algorithm is selected based on the type of atomic orbitals in Wfn.uvsr.
---
## paramters
refWfn::Wfn         -> Wfn object to which the RCCD equations will be applied.

maxit::Int          -> maximum number of coupled cluster iterations.

doprint::Bool=false -> whether or not to print energy and timing information to 
    stdout.
## output
ccenergy::Float -> final RCCD energy. 
"""
function do_rccd(refWfn::Wfn; maxit=40, doprint=false, return_T2=false)
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    epsa = refWfn.epsa
    T = eltype(refWfn.uvsr)
    oovv,ovov,ovvo,oooo,vvvv = make_rccd_integrals(refWfn.uvsr,refWfn.Cao,refWfn.Cav) 
    oovv = CuArray(oovv)
    ovov = CuArray(ovov)
    ovvo = CuArray(ovvo)
    oooo = CuArray(oooo)
    vvvv = CuArray(vvvv)
    T2 = cuzeros(T, nocc, nocc, nvir, nvir)
    T2new = cuzeros(T,nocc,nocc,nvir,nvir)
    Dijab = CuArray(form_Dijab(T2, epsa))
    T2_init!(T2, oovv, Dijab)
    if doprint
        println("@RMP2 ", ccenergy(T2, oovv))
    end
    Fae = form_Fae(T2, oovv)
    Fmi = form_Fmi(T2, oovv)
    Wmnij = form_Wmnij(oooo, oovv, T2)
    Wabef = form_Wabef(vvvv, oovv, T2)
    WmBeJ = form_WmBeJ(ovvo, oovv, T2)
    WmBEj = form_WmBEj(oovv, ovov, T2)
    dt = @elapsed for i in 0:maxit-1 #TODO: implement RMS check
        cciter(
            T2new
            T2,
            oovv,
            vvvv,
            oooo,
            ovov,
            ovvo,
            Dijab,
            Fae,
            Fmi,
            Wabef,
            Wmnij,
            WmBeJ,
            WmBEj,
        )
        if doprint
            println("$i @CCD ", ccenergy(T2, oovv))
        end
        T2 .= T2new
    end
    if doprint
        println("CCD iterations computed in $dt s")
    end
    if return_T2
        return ccenergy(T2, oovv), T2
    else
        return ccenergy(T2, oovv)
    end
end

function make_rccd_integrals(gao::Array,Cao,Cav)
    oovv = tei_transform(gao, Cao, Cav, Cao, Cav, "oovv")
    vvvv = tei_transform(gao, Cav, Cav, Cav, Cav, "vvvv")
    ovvo = tei_transform(gao, Cao, Cav, Cav, Cao, "ovvo")
    ovov = tei_transform(gao, Cao, Cao, Cav, Cav, "oovv")
    oooo = tei_transform(gao, Cao, Cao, Cao, Cao, "oooo")
    oovv = permutedims(oovv,[1,3,2,4])
    ovov = permutedims(ovov,[1,3,2,4])
    ovvo = permutedims(ovvo,[1,3,2,4])
    oooo = permutedims(oooo,[1,3,2,4])
    vvvv = permutedims(vvvv,[1,3,2,4])
    return oovv,ovov,ovvo,oooo,vvvv
end

function ccenergy(tiJaB, ijab)
    #ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    @tensoropt ecc[] := ijab[i,j,a,b]*2*tiJaB[i,j,a,b] - ijab[i,j,a,b]*tiJaB[j,i,a,b]
    return ecc[]
end

@inbounds @fastmath function cciter(
    tiJaB_d,
    tiJaB_i,
    oovv,
    vvvv,
    oooo,
    ovov,
    ovvo,
    Dijab,
    Fae,
    Fmi,
    Wabef,
    Wmnij,
    WmBeJ,
    WmBEj,
)
    form_Wabef!(Wabef, vvvv, oovv, tiJaB_i)
    form_WmBeJ!(WmBeJ, ovvo, oovv, tiJaB_i)
    form_WmBEj!(WmBEj, oovv, ovov, tiJaB_i)
    form_Wmnij!(Wmnij, oooo, oovv, tiJaB_i)
    form_Fae!(Fae, tiJaB_i, oovv)
    form_Fmi!(Fmi, tiJaB_i, oovv)

    tiJaB_d = form_T2(tiJaB_d,tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, oovv, Dijab)
    return tiJaB_d
end

function T2_init!(tiJaB, iajb, Dijab)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    tiJaB .= iajb ./ Dijab
end

function form_Fae(tiJaB, menf)
    dt = eltype(menf)
    nvir = size(tiJaB, 4)
    Fae = cuzeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, menf)
    return Fae
end
function form_Fae!(Fae, tiJaB, mnef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    @tensoropt begin
        Fae[a,e] = -1*(mnef[m,n,e,f]*(2*tiJaB[m,n,a,f] - tiJaB[n,m,a,f]))
    end
end

function form_Fmi(tiJaB, menf)
    dt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Fmi = cuzeros(dt, nocc, nocc)
    form_Fmi!(Fmi, tiJaB, menf)
    return Fmi
end
function form_Fmi!(Fmi, tiJaB, mnef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    @tensoropt begin
        Fmi[m,i] = mnef[m,n,e,f]*(2*tiJaB[i,n,e,f] - tiJaB[i,n,f,e])
    end
    return Fmi
end

function form_T2(tiJaB_d,tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, ijab, Dijab)
    dtt = eltype(tiJaB_i)
    nocc = size(Wmnij, 1)
    nvir = size(tiJaB_i, 4)
    rocc = 1:nocc
    rvir = 1:nvir 
    #tiJaB_d = cuzeros(eltype(tiJaB_i),size(ijab[:,:,:,:]))
    #tiJaB_d = CuArray(tiJaB_d)
    @tensoropt begin
        tiJaB_d[i,j,a,b] = (ijab[i,j,a,b] + tiJaB_i[i,j,a,e]*Fae[b,e] 
                            + tiJaB_i[j,i,b,e]*Fae[a,e]
                           - tiJaB_i[i,m,a,b]*Fmi[m,j] 
                           - tiJaB_i[m,j,a,b]*Fmi[m,i]
                           + tiJaB_i[m,n,a,b]*Wmnij[m,n,i,j]
                           + tiJaB_i[i,j,e,f]*Wabef[a,b,e,f] 
                           + tiJaB_i[i,m,a,e]*WmBeJ[m,b,e,j]*2
                           - tiJaB_i[m,i,a,e]*WmBeJ[m,b,e,j]
                           + tiJaB_i[i,m,a,e]*WmBEj[m,b,e,j]
                           + tiJaB_i[m,i,b,e]*WmBEj[m,a,e,j]
                           + tiJaB_i[m,j,a,e]*WmBEj[m,b,e,i]
                           + tiJaB_i[j,m,b,e]*WmBeJ[m,a,e,i]*2
                           - tiJaB_i[m,j,b,e]*WmBeJ[m,a,e,i]
                           + tiJaB_i[j,m,b,e]*WmBEj[m,a,e,i])

    end
    tiJaB_d .= tiJaB_d ./ Dijab
    return tiJaB_d
end
function form_Wmnij(minj, menf, tiJaB)
    dtt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Wmnij = cuzeros(dtt, nocc, nocc, nocc, nocc)
    form_Wmnij!(Wmnij, minj, menf, tiJaB)
    return Wmnij
end
function form_Wmnij!(Wmnij, mnij, mnef, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    @tensoropt begin
        Wmnij[m,n,i,j] = mnij[m,n,i,j] + tiJaB[i,j,e,f]*mnef[m,n,e,f]/2
    end
end

function form_Wabef(aebf, mnef, tiJaB)
    dt = eltype(aebf)
    nvir = size(tiJaB, 4)
    Wabef = cuzeros(dt, nvir, nvir, nvir, nvir)
    form_Wabef!(Wabef, aebf, mnef, tiJaB)
    return Wabef
end
function form_Wabef!(Wabef, abef, mnef, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    @tensoropt begin
        Wabef[a,b,e,f] = tiJaB[m,n,a,b]*mnef[m,n,e,f]/2 + abef[a,b,e,f]
    end
end

function form_WmBeJ(mebj, iajb, tiJaB)
    dtt = eltype(iajb)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBeJ = cuzeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBeJ!(WmBeJ, mebj, iajb, tiJaB)
    return WmBeJ
end
function form_WmBeJ!(WmBeJ, mbej, mnef, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    @tensoropt begin
        WmBeJ[m,b,e,j] = (mbej[m,b,e,j] + mnef[m,n,e,f]*(2*tiJaB[n,j,f,b] - tiJaB[j,n,f,b])/2
                          - mnef[n,m,e,f]*tiJaB[n,j,f,b]/2)
    end
end


function form_WmBEj(nemf, mjbe, tiJaB)
    dtt = eltype(nemf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBEj = cuzeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBEj!(WmBEj, nemf, mjbe, tiJaB)
    return WmBEj
end
function form_WmBEj!(WmBEj, nmef, mbje, tiJaB)
    dtt = eltype(WmBEj)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    @tensoropt begin
        WmBEj[m,b,e,j] = -mbje[m,b,j,e] + tiJaB[j,n,f,b]*nmef[n,m,e,f]/2.0
    end
    return WmBEj
end
end #module
