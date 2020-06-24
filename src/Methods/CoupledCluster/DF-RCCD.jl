module DFRCCD
using Fermi.DF
using Fermi.Wavefunction
using TensorOperations
using Fermi.Transformation
include("Denominators.jl")
export do_df_rccd
"""
    do_df_rccd
"""
function do_df_rccd(refWfn::Wfn; maxit=40, doprint=false, return_T2=false, dfbname="default")
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    epsa = refWfn.epsa
    T = eltype(refWfn.uvsr)
    bov,bvo,boo,bvv = make_df_rccd_integrals(refWfn::Wfn; dfbname=dfbname)
    T2 = zeros(T, nocc, nocc, nvir, nvir)
    Dijab = form_Dijab(T2, epsa)
    T2_init!(T2, bov, Dijab)
    ecc = ccenergy(T2, bov)
    if doprint println("@DF-RMP2 $ecc") end
    WmBeJ = form_WmBeJ(bov, bvo, T2)
    WmBEj = form_WmBEj(boo, bvv, bov, T2)
    Wmnij = form_Wmnij(boo, bov, T2)
    Fae = form_Fae(T2, bov)
    Fmi = form_Fmi(T2, bov)
    dt = @elapsed for i in 0:maxit-1 #TODO: implement RMS check
        T2 = cciter(
            T2,
            bov,
            bvo,
            boo,
            bvv,
            Dijab,
            Fae,
            Fmi,
            Wmnij,
            WmBeJ,
            WmBEj,
        )
        if doprint
            println("$i @CCD ", ccenergy(T2, bov))
        end
    end
    if doprint
        println("CCD iterations computed in $dt s")
    end
    if return_T2
        return ccenergy(T2, bov), T2
    else
        return ccenergy(T2, bov)
    end
    return ecc
end
@inbounds @fastmath function cciter(
    tiJaB_i,
    bov,
    bvo,
    boo,
    bvv,
    Dijab,
    Fae,
    Fmi,
    Wmnij,
    WmBeJ,
    WmBEj,
)
    form_WmBeJ!(WmBeJ, bov, bvo, tiJaB_i)
    form_WmBEj!(WmBEj, boo, bvv, bov, tiJaB_i)
    form_Wmnij!(Wmnij, boo, bov, tiJaB_i)
    form_Fae!(Fae, tiJaB_i, bov)
    form_Fmi!(Fmi, tiJaB_i, bov)

    tiJaB_d = form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wmnij, bov, bvv, Dijab)
    return tiJaB_d
end

function make_df_rccd_integrals(refWfn::Wfn; dfbname="default")
    bμν = DF.make_bμν(refWfn; dfbname=dfbname)
    Cao = refWfn.Cao
    Cav = refWfn.Cav
    @tensoropt begin
        bov[i,a,Q] := Cao[μ,i]*bμν[μ,ν,Q]*Cav[ν,a]
        bvo[a,i,Q] := Cav[μ,a]*bμν[μ,ν,Q]*Cao[ν,i]
        boo[i,j,Q] := Cao[μ,i]*bμν[μ,ν,Q]*Cao[ν,j]
        bvv[a,b,Q] := Cav[μ,a]*bμν[μ,ν,Q]*Cav[ν,b]
    end
    return bov,bvo,boo,bvv
end
function ccenergy(tiJaB, bov)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    for i in rocc
        for j in rocc
            bi = bov[i,:,:]
            bj = bov[j,:,:]
            @tensor begin
                bAB[a,b] := bi[a,Q]*bj[b,Q]
            end
            for a in rvir
                for b in rvir
                    ecc += bAB[a,b] * 2 * tiJaB[i,j,a,b]
                    ecc -= bAB[a,b] * tiJaB[j,i,a,b]
                end
            end
        end
    end
    return ecc
end
function T2_init!(tiJaB,bov,Dijab)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    nao  = size(bov,   3)
    rocc = 1:nocc
    rvir = 1:nvir
    #@tensor begin
    #    tiJaB[i,j,a,b] = bov[i,a,Q]*bov[j,b,Q]#/Dijab[i,j,a,b]
    #end
    for i in rocc
        for j in rocc
            bi = bov[i,:,:]
            bj = bov[j,:,:]
            @tensor begin
                bAB[a,b] := bi[a,Q]*bj[b,Q]
            end
            for a in rvir
                for b in rvir
                    tiJaB[i,j,a,b] = bAB[a,b]/Dijab[i,j,a,b]
                end
            end
        end
    end
end
function form_Fae(tiJaB, bov)
    dt = eltype(bov)
    nvir = size(tiJaB, 4)
    Fae = zeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, bov)
    return Fae
end
function form_Fae!(Fae,tiJaB, bov)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    @tensoropt begin
        Fae[a,e] = -1*bov[m,e,Q]*bov[n,f,Q]*(2*tiJaB[m,n,a,f] - tiJaB[n,m,a,f])
    end
    return Fae
end
function form_Fmi(tiJaB, bov)
    dt = eltype(bov)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Fmi = zeros(dt, nocc, nocc)
    form_Fmi!(Fmi, tiJaB, bov)
    return Fmi
end
function form_Fmi!(Fmi, tiJaB, bov)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    @tensoropt begin
        Fmi[m,i] = bov[m,e,Q]*bov[n,f,Q]*(2*tiJaB[i,n,e,f] - tiJaB[i,n,f,e])
    end
    return Fmi
end
function form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wmnij, bov, bvv, Dijab)
    dtt = eltype(tiJaB_i)
    nocc = size(Wmnij, 1)
    nvir = size(tiJaB_i, 4)
    rocc = 1:nocc
    rvir = 1:nvir 
    #ijab = permutedims(iajb,[1,3,2,4])
    tiJaB_d = zeros(size(tiJaB_i[:,:,:,:]))
    #for i in rocc
    #    for j in rocc
    #        _T2 = tiJaB_i[i,j,:,:]
    #for a in rvir
    #    for b in rvir
    #        bavQ = bvv[a,:,:]
    #        bbvQ = bvv[b,:,:]
    #        for i in rocc
    #            for j in rocc
    #                tiJaB_i_ij = tiJaB_i[i,j,:,:]
    #                @tensor begin
    #                    temp[] := tiJaB_i_ij[e,f]*bavQ[e,Q]*bbvQ[f,Q]
    #                end
    #                tiJaB_d[i,j,a,b] = temp[]
    #            end
    #        end
    #    end
    #end
    @tensoropt begin
        tiJaB_d[i,j,a,b] += (bov[i,a,Q]*bov[j,b,Q] + tiJaB_i[i,j,a,e]*Fae[b,e] + tiJaB_i[j,i,b,e]*Fae[a,e]
                           - tiJaB_i[i,m,a,b]*Fmi[m,j] - tiJaB_i[m,j,a,b]*Fmi[m,i]
                           + tiJaB_i[m,n,a,b]*Wmnij[m,n,i,j]
                           #+ tiJaB_i[i,j,e,f]*Wabef[a,b,e,f] 
                           + tiJaB_i[i,j,e,f]*(
                                               tiJaB_i[m,n,a,b]
                                               *bov[m,e,Q]*bov[n,f,Q]/2)
                                               #+bvv[a,e,Q]*bvv[b,f,Q])
                           + tiJaB_i[i,j,e,f]*bvv[a,e,Q]*bvv[b,f,Q]
                           + tiJaB_i[i,m,a,e]*WmBeJ[m,b,e,j]*2
                           - tiJaB_i[m,i,a,e]*WmBeJ[m,b,e,j]
                           + tiJaB_i[i,m,a,e]*WmBEj[m,b,e,j]
                           + tiJaB_i[m,i,b,e]*WmBEj[m,a,e,j]
                           + tiJaB_i[m,j,a,e]*WmBEj[m,b,e,i]
                           + tiJaB_i[j,m,b,e]*WmBeJ[m,a,e,i]*2
                           - tiJaB_i[m,j,b,e]*WmBeJ[m,a,e,i]
                           + tiJaB_i[j,m,b,e]*WmBEj[m,a,e,i])

    tiJaB_d = tiJaB_d ./ Dijab
    end
end
function form_Wmnij(boo, bov, tiJaB)
    dtt = eltype(bov)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Wmnij = zeros(dtt, nocc, nocc, nocc, nocc)
    form_Wmnij!(Wmnij, boo, bov, tiJaB)
    return Wmnij
end
function form_Wmnij!(Wmnij, boo, bov, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    @tensoropt begin
        Wmnij[m,n,i,j] = boo[m,i,Q]*boo[n,j,Q] + tiJaB[i,j,e,f]*bov[m,e,Q]*bov[n,f,Q]/2
    end
end
function form_WmBeJ(bov, bvo, tiJaB)
    dtt = eltype(bov)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBeJ = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBeJ!(WmBeJ, bov, bvo, tiJaB)
    return WmBeJ
end
function form_WmBeJ!(WmBeJ, bov, bvo, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    @tensoropt begin
        WmBeJ[m,b,e,j] = (bov[m,e,Q]*bvo[b,j,Q] + bov[m,e,Q]*bov[n,f,Q]*(
                            2*tiJaB[n,j,f,b] - tiJaB[j,n,f,b])/2
                          - bov[m,e,Q]*bov[n,f,Q]*tiJaB[n,j,f,b]/2)
    end
end
function form_WmBEj(boo,bvv,bov,tiJaB)
    dtt = eltype(boo)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBEj = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBEj!(WmBEj, boo, bvv, bov, tiJaB)
    return WmBEj
end
function form_WmBEj!(WmBEj, boo, bvv, bov, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    @tensoropt begin
        WmBEj[m,b,e,j] = ( -boo[m,j,Q]*bvv[b,e,Q] 
                          + tiJaB[j,n,f,b]*bov[n,e,Q]*bov[m,f,Q]/2)
    end
end
end #module
