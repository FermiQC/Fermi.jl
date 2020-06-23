module ROCCD

using Fermi.Wavefunction
using Fermi.Transformation
using TensorOperations
using LinearAlgebra

include("Denominators.jl")

export do_roccd

function do_roccd(refWfn::Wfn,maxit; doprint::Bool=false, return_T2::Bool=false)
    nocca = refWfn.nalpha
    noccb = refWfn.nbeta
    nvira = refWfn.nvira
    nvirb = refWfn.nvirb
    N = max(nocca,noccb)
    Δ = nocca - noccb
    uvsr = refWfn.uvsr
    Ca = refWfn.Ca
    Cb = refWfn.Cb
    Cao = refWfn.Cao
    Cbo = refWfn.Cbo
    Cav = refWfn.Cav
    Cbv = refWfn.Cbv
    h = refWfn.hao
    @tensor begin
        ha[p,q] := Ca[m,q]*h[m,n]*Ca[n,p]
    end
    @tensor begin
        hb[p,q] := Cb[m,q]*h[m,n]*Ca[n,p]
    end
    xoxo = permutedims(tei_transform(uvsr,Ca,Ca,Cao,Cao,"xoxo"),[1,3,2,4])
    xxoo = permutedims(tei_transform(uvsr,Ca,Cao,Ca,Cao,"xxoo"),[1,3,2,4])
    XOXO = permutedims(tei_transform(uvsr,Cb,Cb,Cbo,Cbo,"XOXO"),[1,3,2,4])
    fa = form_fa(ha,xoxo,xxoo,XOXO)
    xxoo = nothing
    XXOO = permutedims(tei_transform(uvsr,Cb,Cbo,Cb,Cbo,"XXOO"),[1,3,2,4])
    fb = form_fb(hb,XOXO,XXOO,xoxo)
    xoxo = nothing
    XOXO = nothing
    XXOO = nothing

    oOvV = permutedims(tei_transform(uvsr,Cao,Cav,Cbo,Cbv,"oOvV"),[1,3,2,4])
    oovv = permutedims(tei_transform(uvsr,Cao,Cav,Cao,Cav,"oovv"),[1,3,2,4])
    oovv = oovv - permutedims(oovv,[1,2,4,3])
    OOVV = permutedims(tei_transform(uvsr,Cbo,Cbv,Cbo,Cbv,"OOVV"),[1,3,2,4])
    OOVV = OOVV - permutedims(OOVV,[1,2,4,3])

    oOoO = permutedims(tei_transform(uvsr,Cao,Cao,Cbo,Cbo,"oOoO"),[1,3,2,4])
    oooo = permutedims(tei_transform(uvsr,Cao,Cao,Cao,Cao,"oooo"),[1,3,2,4])
    oooo = oooo - permutedims(oooo,[1,2,4,3])
    OOOO = permutedims(tei_transform(uvsr,Cbo,Cbo,Cbo,Cbo,"OOOO"),[1,3,2,4])
    OOOO = OOOO - permutedims(OOOO,[1,2,4,3])
    vVvV = permutedims(tei_transform(uvsr,Cav,Cav,Cbv,Cbv,"vVvV"),[1,3,2,4])
    vvvv = permutedims(tei_transform(uvsr,Cav,Cav,Cav,Cav,"vvvv"),[1,3,2,4])
    oVvO = permutedims(tei_transform(uvsr,Cao,Cav,Cbv,Cbo,"oVvO"),[1,3,2,4])
    vOvO = permutedims(tei_transform(uvsr,Cav,Cav,Cbo,Cbo,"vOvO"),[1,3,2,4])
    oVoV = permutedims(tei_transform(uvsr,Cao,Cao,Cbv,Cbv,"oVoV"),[1,3,2,4])
    vOoV = permutedims(tei_transform(uvsr,Cav,Cao,Cbo,Cbv,"vOoV"),[1,3,2,4])

    Dijab = form_Dijab(oovv,fa)
    DiJaB = form_DiJaB(oOvV,fa,fb)
    DIJAB = form_Dijab(OOVV,fb)

    tijab = zeros(size(oovv))
    tiJaB = zeros(size(oOvV))
    tIJAB = zeros(size(OOVV))

    tijab .= oovv ./ Dijab
    tiJaB .= oOvV ./ DiJaB
    tIJAB .= OOVV ./ DIJAB

    ecc = ccenergy(oOvV,oovv,OOVV,tijab,tiJaB,tIJAB)
    println("@RO-MP2 $ecc")

    Fae = zeros(nvira,nvira)
    FAE = zeros(nvirb,nvirb)
    Fmi = zeros(nocca,nocca)
    FMI = zeros(noccb,noccb)
    Fme = fa[1:nocca,nocca+1:nocca+nvira]
    FME = fb[1:noccb,noccb+1:noccb+nvirb]
    Wmnij = zeros(nocca,nocca,nocca,nocca)
    WmNiJ = zeros(nocca,noccb,nocca,noccb)
    WMniJ = zeros(noccb,nocca,nocca,noccb)
    WMNIJ = zeros(noccb,noccb,noccb,noccb)
    Wabef = zeros(nvira,nvira,nvira,nvira)
    WaBeF = zeros(nvira,nvirb,nvira,nvirb)
    WABEF = zeros(nvirb,nvirb,nvirb,nvirb)
    Wmbej = zeros(nocca,nvira,nvira,nocca)
    WmBeJ = zeros(nocca,nvirb,nvira,noccb)
    WmBEj = zeros(nocca,nvirb,nvirb,nocca)
    WMBEJ = zeros(noccb,nvirb,nvirb,noccb)
    WMbEj = zeros(noccb,nvira,nvirb,nocca)
    WMbeJ = zeros(noccb,nvira,nvira,noccb)


    form_Fae!(Fae,fa[nocca+1:nocca+nvira,nocca+1:nocca+nvira],oOvV,oovv,tijab,tiJaB)
    form_FAE!(FAE,fb[noccb+1:nocca+nvira,noccb+1:noccb+nvirb],oOvV,OOVV,tIJAB,tiJaB)
    form_Fmi!(Fmi,fa[1:nocca,1:nocca],oovv,oOvV,tijab,tiJaB)
    form_FMI!(FMI,fb[1:noccb,1:noccb],OOVV,oOvV,tIJAB,tiJaB)
    form_Wmnij!(Wmnij,oooo,oovv,tijab)
    form_WmNiJ!(WmNiJ,oOoO,oOvV,tiJaB)
    form_WMNIJ!(WMNIJ,OOOO,OOVV,tIJAB)
    form_WMniJ!(WMniJ,oOoO,oOvV,tiJaB)
    form_Wabef!(Wabef,vVvV,oOvV,tijab)
    form_WaBeF!(WaBeF,vVvV,oOvV,tiJaB)
    form_WABEF!(WABEF,vVvV,oOvV,tIJAB)
    form_Wmbej!(Wmbej,oVvO,vOvO,oOvV,tijab,tiJaB)
    form_WmBeJ!(WmBeJ,oVvO,oOvV,tIJAB,tiJaB)
    form_WmBEj!(WmBEj,oVoV,oOvV,tiJaB)
    form_WMBEJ!(WMBEJ,oVvO,vOvO,oOvV,tIJAB,tiJaB)
    form_WMbEj!(WMbEj,vOoV,oOvV,tijab,tiJaB)
    form_WMbeJ!(WMbeJ,vOvO,oOvV,tiJaB)
end
function ccenergy(oOvV,oovv,OOVV,tijab,tiJaB,tIJAB)
    @tensoropt begin
        ecc1[] := ( (1/4)*tijab[i,j,a,b]*oovv[i,j,a,b])
           #       + tiJaB[i,j,a,b]*oOvV[i,j,a,b]
           #       + (1/4)*tIJAB[i,j,a,b]*(oOvV[i,j,a,b] - oOvV[j,i,a,b]))
    end
    @tensoropt begin
        ecc2[] := (tiJaB[i,j,a,b]*oOvV[i,j,a,b])
    end
    @tensoropt begin
        ecc3[] := ((1/4)*tIJAB[i,j,a,b]*OOVV[i,j,a,b])
    end
    ecc1 = ecc1[]
    ecc2 = ecc2[]
    ecc3 = ecc3[]
    println("AA $ecc1")
    println("AB $ecc2")
    println("BB $ecc3")
    return ecc1 + ecc2 + ecc3
end
function form_tijab(oovv,Fae,Fmi,Wmnij,Wabef,Wmbej,WMbEj,tijab,tiJaB,Dijab)
    @tensoropt begin
        _tijab[i,j,a,b] := (oovv[i,j,a,b]
                            + tijab[i,j,a,e]*Fae[b,e]
                            - tijab[i,j,b,e]*Fae[a,e]
                            - tijab[i,m,a,b]*Fmi[m,j]
                            + tijab[j,m,a,b]*Fmi[m,i]
                            + (1/2)*tijab[m,n,a,b]*Wmnij[m,n,i,j]
                            + (1/2)*tijab[i,j,e,f]*Wabef[a,b,e,f]
                            + tijab[i,m,a,e]*Wmbej[m,b,e,j]
                            + tiJaB[i,m,a,e]*WMbEj[m,b,e,j]
                            - tijab[i,m,b,e]*Wmbej[m,a,e,j]
                            - tiJaB[i,m,b,e]*WMbEj[m,a,e,j]
                            - tijab[j,m,a,e]*Wmbej[m,b,e,i]
                            - tiJaB[j,m,a,e]*WMbEj[m,b,e,i]
                            + tijab[j,m,b,e]*Wmbej[m,a,e,i]
                            + tiJaB[j,m,b,e]*WMbEj[m,a,e,i])
    end
    return _tijab ./ Dijab
end
function form_tiJaB(oOvV,Fae,FAE,Fmi,FMI,WmNiJ,WaBeF,WMBEJ,WMbEj,Wmbej,WmBEj,WMbeJ,tiJab,tijab,tIJAB)
    @tensoropt begin
        _tiJaB[i,j,a,b] := (oOvV[i,j,a,b] 
                            + tiJaB[i,j,a,e]*FAE[b,e]
                            + tiJaB[i,j,e,b]*Fae[a,e]
                            - tiJaB[i,m,a,b]*FMI[m,j]
                            - tiJaB[m,j,a,b]*Fmi[m,i]
                            + (1/2)*tiJaB[m,n,a,b]*WmNiJ[m,n,i,j]
                            - (1/2)*tiJaB[n,m,a,b]*WmNiJ[m,n,i,j] #?
                            + (1/2)*tiJaB[i,j,e,f]*WaBeF[a,b,e,f]
                            - (1/2)*tiJaB[i,j,f,e]*WaBeF[a,b,e,f]
                            + tijab[i,m,a,e]*WmBeJ[m,b,e,j]
                            + tiJaB[i,m,a,e]*WMBEJ[m,b,e,j]
                            + tiJaB[i,m,e,b]*WMbeJ[m,a,e,j]
                            + tiJaB[m,j,a,e]*WmBEj[m,b,e,i]
                            + tIJAB[j,m,b,e]*WMbEj[m,a,e,i]
                            + tiJaB[m,j,e,b]*Wmbej[m,a,e,i])
    end
    return _tiJaB ./ Dijab
end
#function form_tIJAB(oOvV,FAE,FMI,WMNIJ,WABEF,WMBEJ,WmBeJ,tiJaB,tIJAB,DIJAB)
#    @tensoropt begin
#        _tIJAB[i,j,a,b] := (oOvV[i,j,a,b] - oOvV[j,i,a,b]
#                            + tIJAB[i,j,a,e]*FAE[b,e]
#                            - tIJAB[i,j,b,e]*FAE[a,e]
#                            - tIJAB[i,m,a,b]*FMI[m,j]
#                            + tIJAB[j,m,a,b]*FMI[m,i]
#                            + (1/2)*tIJAB[m,n,a,b]*WMNIJ[m,n,i,j]
#                            + (1/2)*tIJAB[i,j,e,f]*WABEF[a,b,e,f]
#                            + tIJAB[i,m,a,e]*WMBEJ[m,b,e,j]
#                            + tiJaB[m,i,e,a]*WmBeJ[m,b,e,j]
#                            - tIJAB[i,m,b,e]*WMBEJ[m,a,e,j]
#                            - tiJaB[m,i,e,b]*WmBeJ[m,a,e,j]
#                            - tIJAB[j,m,a,e]*WMBEJ[m,b,e,i]
#                            - tiJaB[m,j,e,a]*WmBeJ[m,b,e,i]
#                            + tIJAB[j,m,b,e]*WMBEJ[m,a,e,i]
#                            = tiJaB[m,j,e,b]*WmBeJ[m,a,e,i])
#    end
#    return _tIJAB ./ DIJAB
#end
function form_fa(ha,xoxo,xxoo,XOXO)
    fa = zeros(size(ha))
    @tensor begin
        fa[p,q] = ha[p,q] + xoxo[p,k,q,k] - xxoo[p,q,k,k] + XOXO[p,K,q,K]
    end
end
function form_fb(hb,XOXO,XXOO,xoxo)
    fb = zeros(size(hb))
    @tensor begin
        fb[p,q] = hb[p,q] + XOXO[p,k,q,k] - XXOO[p,q,k,k] + xoxo[p,K,q,K]
    end
end
function form_Fae!(Fae,fa,oOvV,oovv,tijab,tiJaB)
    @tensoropt begin
        Fae[a,e] := (-0.5*tijab[m,n,a,f]*oovv[m,n,e,f])
        Fae[a,e] -= (-1*tiJaB[m,n,a,f]*oOvV[m,n,e,f])
        #Fae[a,e] -= (-0.5*tIJAB[m,n,a,f]*OOVV[m,n,e,f])
    end
    Fae += fa - diagm(diag(fa))
    return Fae
end
function form_FAE!(FAE,fb,oOvV,OOVV,tIJAB,tiJaB)
    @tensoropt begin
        Fae[a,e] := (-0.5*tIJAB[m,n,a,f]*OOVV[m,n,e,f])
        Fae[a,e] -= (-1*tiJaB[n,m,f,a]*oOvV[n,m,f,e])
        Fae[a,e] -= (-0.5*tIJAB[m,n,a,f]*OOVV[m,n,e,f])
    end
    FAE += fb - diagm(diag(fb))
    return FAE
end
function form_Fmi!(Fmi,fa,oovv,oOvV,tijab,tiJaB)
    @tensoropt begin
        Fmi[m,i] = (0.5*(tijab[i,n,e,f]*(oovv[m,n,e,f]))
                     + tiJaB[i,n,e,f]*oOvV[m,n,e,f])
    end
    Fmi += fa - diagm(diag(fa))
    return Fmi
end
function form_FMI!(FMI,fb,OOVV,oOvV,tIJAB,tiJaB)
    @tensoropt begin
        FMI[m,i] = (0.5*tIJAB[i,n,e,f]*OOVV[m,n,e,f]
                    + tiJaB[n,i,f,e]*oOvV[n,m,f,e])
    end
    FMI += fb - diagm(diag(fb))
    return FMI
end
function form_Wmnij!(Wmnij,oooo,oovv,tijab)
    @tensoropt begin
        Wmnij[m,n,i,j] = (oooo[m,n,i,j]
                          + (1/4)*(tijab[i,j,e,f]
                                   *(oovv[m,n,e,f])))
    end
    return Wmnij
end
function form_WMNIJ!(WMNIJ,OOOO,OOVV,tIJAB)
    @tensoropt begin
        WMNIJ[m,n,i,j] = (OOOO[m,n,i,j]
                          + (1/4)*tIJAB[i,j,e,f]*OOVV[m,n,e,f])
    end
    return WMNIJ
end
function form_WmNiJ!(WmNiJ,oOoO,oOvV,tiJaB)
    @tensoropt begin
        WmNiJ[m,n,i,j] = (oOoO[m,n,i,j] + (1/4)*tiJaB[i,j,e,f]*oOvV[m,n,e,f]
                          + (1/4)*tiJaB[i,j,f,e]*oOvV[m,n,f,e])
    end
    return WmNiJ
end
function form_WMniJ!(WMniJ,oOoO,oOvV,tiJaB)
    println(size(tiJaB))
    println(size(oOvV))
    @tensoropt begin
        WMniJ[m,n,i,j] = (-1*oOoO[n,m,i,j]
                          - (1/4)*tiJaB[i,j,e,f]*oOvV[m,n,f,e])
                          #- (1/4)*tiJaB[i,j,f,e]*oOvV[m,n,e,f])
    end
    return WMniJ
end
function form_Wabef!(Wabef,vVvV,oOvV,tijab)
    @tensoropt begin
        Wabef[a,b,e,f] = (vVvV[a,b,e,f] - vVvV[b,a,e,f] 
                          + (1/4)*tijab[m,n,a,b]
                          *(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
    end
end
function form_WaBeF!(WaBeF,vVvV,oOvV,tiJaB)
    @tensoropt begin
        WaBeF[a,b,e,f] = vVvV[a,b,e,f] + (1/2)*tiJaB[m,n,a,b]*oOvV[m,n,e,f]
    end
end
function form_Wmbej!(Wmbej,oVvO,vOvO,oOvV,tijab,tiJaB)
    @tensoropt begin
        Wmbej[m,b,ej] = (oVvO[m,b,e,j] - vOvO[b,m,e,j] 
                         - (1/2)*tijab[j,n,f,b]
                         *(oOvV[m,n,e,f] - oOvV[n,m,e,f])
                         + (1/2)*tiJaB[j,n,b,f]*oOvV[m,n,e,f])
    end
    return Wmbej
end
function form_WmBeJ!(WmBeJ,oVvO,oOvV,tIJAB,tiJaB)
    @tensoropt begin
        WmBeJ[m,b,e,j] = (oVvO[m,b,e,j] 
                          - (1/2)*tIJAB[j,n,f,b]*oOvV[m,n,e,f]
                          + (1/2)*tiJaB[n,j,f,b]
                            *(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
    end
    return WmBeJ
end
function form_WmBEj!(WmBEj,oVoV,oOvV,tiJaB)
    @tensoropt begin
        WmBEj[m,b,e,j] = (-1*oVoV[m,b,j,e]
                          + (1/2)*tiJaB[j,n,f,b]*oOvV[m,n,f,e])
    end
    return WmBEj
end
function form_WMBEJ!(WMBEJ,oVvO,vOvO,oOvV,tIJAB,tiJaB)
    @tensoropt begin
        WMBEJ[m,b,e,j] = (oVvO[m,b,e,j] - vOvO[b,m,e,j]
                          - (1/2)*tIJAB[j,n,f,b]
                            *(oOvV[m,n,e,f] - oOvV[n,m,e,f])
                            + (1/2)*tiJaB[n,j,f,b]*oOvV[n,m,f,e])
    end
    return WMBEJ
end
function form_WMbEj!(WMbEj,vOoV,oOvV,tijab,tiJaB)
    @tensoropt begin
        WMbEj[m,b,e,j] = (vOoV[b,m,j,e]
                          - (1/2)*tijab[j,n,f,b]*oOvV[n,m,f,e]
                          + (1/2)*tiJaB[j,n,b,f]
                            *(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
    end
    return WMbEj
end
function form_WMbeJ(WMbeJ,vOvO,oOvV,tiJaB)
    @tensoropt begin
        WMbeJ[m,b,e,j] = (-1*vOvO[b,m,e,j]
                          + (1/2)*tiJaB[n,j,b,f]*oOvV[n,m,e,f])
    end
    return WMbeJ
end
end #module
