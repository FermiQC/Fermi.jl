module UCCSD

using Fermi.Wavefunction
using Fermi.Transformation
using TensorOperations
using LinearAlgebra

include("Denominators.jl")

export do_uccsd

function do_uccsd(refWfn::Wfn; maxit=40, doprint::Bool=false, return_T2::Bool=false)
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
        hb[p,q] := Cb[m,q]*h[m,n]*Cb[n,p]
    end


    oooo = tei_transform(uvsr,Cao,Cao,Cao,Cao,"oooo")
    OOOO = tei_transform(uvsr,Cbo,Cbo,Cbo,Cbo,"OOOO")
    ooOO = tei_transform(uvsr,Cao,Cao,Cbo,Cbo,"ooOO")

    escf = refWfn.vnuc
    rocca = 1:nocca
    roccb = 1:noccb
    for i in rocca
        escf += ha[i,i]
    end
    for i in roccb
        escf += ha[i,i]
    end
    for i in rocca
        for j in rocca
            escf += 0.5*oooo[i,i,j,j]
            escf -= 0.5*oooo[i,j,j,i]
        end
    end
    for i in roccb
        for j in roccb
            escf += 0.5*OOOO[i,i,j,j]
            escf -= 0.5*OOOO[i,j,j,i]
        end
    end
    for i in rocca
        for j in roccb
            escf += ooOO[i,i,j,j]
        end
    end
    if doprint println("@UHF $escf") end

    xxoo = tei_transform(uvsr,Ca,Ca,Cao,Cao,"xxoo")
    xoxo = tei_transform(uvsr,Ca,Cao,Ca,Cao,"xoxo")
    XXOO = tei_transform(uvsr,Cb,Cb,Cbo,Cbo,"XXOO")
    XOXO = tei_transform(uvsr,Cb,Cbo,Cb,Cbo,"XOXO")
    xxOO = tei_transform(uvsr,Ca,Ca,Cbo,Cbo,"xxOO")
    XXoo = tei_transform(uvsr,Cb,Cb,Cao,Cao,"XXoo")

    #these dont work, so just creating canonical fock matrix
    fa = form_fa(ha,xxoo,xoxo,xxOO)
    fb = form_fb(hb,XXOO,XOXO,XXoo)

    oovv = permutedims(tei_transform(uvsr,Cao,Cav,Cao,Cav,"oovv"),[1,3,2,4])# - permutedims(permutedims(iajb,[1,4,3,2]),[1,3,2,4])
    IAJB = tei_transform(uvsr,Cbo,Cbv,Cbo,Cbv,"IAJB")
    OOVV = permutedims(IAJB,[1,3,2,4])# - permutedims(permutedims(IAJB,[1,4,3,2]),[1,3,2,4])
    oOvV = permutedims(tei_transform(uvsr,Cao,Cav,Cbo,Cbv,"oOvV"),[1,3,2,4])
    OoVv = permutedims(tei_transform(uvsr,Cbo,Cbv,Cao,Cav,"OoVv"),[1,3,2,4])
    ooov = permutedims(tei_transform(uvsr,Cao,Cao,Cao,Cav,"ooov"),[1,3,2,4])
    ooov = ooov# - permutedims(ooov,[1,2,4,3])
    OOOV = permutedims(tei_transform(uvsr,Cbo,Cbo,Cbo,Cbv,"OOOV"),[1,3,2,4])
    OoOv = permutedims(tei_transform(uvsr,Cbo,Cbo,Cao,Cav,"OoVv"),[1,3,2,4])
    OOVO = permutedims(tei_transform(uvsr,Cbo,Cbv,Cbo,Cbo,"OOVO"),[1,3,2,4])
    oOoV = permutedims(tei_transform(uvsr,Cao,Cao,Cbo,Cbv,"oOvV"),[1,3,2,4])
    oovo = permutedims(tei_transform(uvsr,Cao,Cav,Cao,Cao,"oovo"),[1,3,2,4])
    vovv = permutedims(tei_transform(uvsr,Cav,Cav,Cao,Cav,"vovv"),[1,3,2,4]) 
    vovv = vovv# - permutedims(vovv,[1,2,4,3])
    vOvV = permutedims(tei_transform(uvsr,Cav,Cav,Cbo,Cbv,"vOvV"),[1,3,2,4])
    VoVv = permutedims(tei_transform(uvsr,Cbv,Cbv,Cao,Cav,"VoVv"),[1,3,2,4])
    #there may be a bug in how I antisymmetrized some integrals. 
    VOVV = permutedims(tei_transform(uvsr,Cbv,Cbv,Cbo,Cbv,"VOVV"),[1,3,2,4])
    VOVV = VOVV# - permutedims(VOVV,[1,2,4,3])


    Dia   = form_Dia(oovv,fa)
    DIA   = form_Dia(OOVV,fb)
    Dijab = form_Dijab(oovv,fa)
    DIJAB = form_Dijab(OOVV,fb)
    DiJaB = form_DiJaB(oOvV,fa,fb)

    tia = fa[1:nocca,nocca+1:size(fa)[1]] ./ Dia
    tIA = fb[1:noccb,noccb+1:size(fb)[1]] ./ DIA
    tijab = (oovv - permutedims(oovv,[1,2,4,3])) ./ Dijab
    tIJAB = (OOVV - permutedims(OOVV,[1,2,4,3])) ./ DIJAB
    tiJaB = oOvV ./ DiJaB

    fa_oo = fa[1:nocca,1:nocca]
    fb_oo = fb[1:noccb,1:noccb]
    fa_ov = fa[1:nocca,nocca+1:size(fa)[1]]
    fb_ov = fb[1:noccb,noccb+1:size(fb)[1]]
    fa_vv = fa[nocca+1:size(fa)[1],nocca+1:size(fa)[1]]
    fb_vv = fb[noccb+1:size(fb)[1],noccb+1:size(fb)[1]]
    @tensor begin
        tiatia[m,n,a,f] := tia[m,a]*tia[n,f]
        tiatIA[m,n,a,f] := tia[m,a]*tIA[n,f]
        tIAtIA[m,n,a,f] := tIA[m,a]*tIA[n,f]
        tIAtia[m,n,a,f] := tIA[m,a]*tia[n,f]
        emp2_aa[] := (1/4)*tijab[i,j,a,b]*(oovv[i,j,a,b] - oovv[i,j,b,a])
        emp2_bb[] := (1/4)*tIJAB[i,j,a,b]*(OOVV[i,j,a,b] - OOVV[i,j,b,a])
        emp2_ab[] := tiJaB[i,j,a,b]*oOvV[i,j,a,b]
    end
    emp2 = emp2_aa[] + emp2_bb[] + emp2_ab[]
    if doprint println("@UMP2 $emp2") end

    Fae = form_Fae(fa_vv,fa_ov,oovv,oOvV,vovv,vOvV,tia,tIA,tiatia,tiatIA,tijab,tiJaB)
    FAE = form_FAE(fb_vv,fb_ov,OOVV,OoVv,VOVV,VoVv,tia,tIA,tIAtIA,tIAtia,tIJAB,tiJaB)
    Fmi = form_Fmi(fa_oo,fa_ov,ooov,oovo,oOoV,oovv,oOvV,tia,tIA,tiatia,tiatIA,tijab,tiJaB)
    FMI = form_FMI(fb_oo,fb_ov,OOOV,OOVO,OoOv,OOVV,OoVv,tia,tIA,tIAtIA,tIAtia,tIJAB,tiJaB)
    println(FMI)
    return emp2

end
function form_fa(ha,xxoo,xoxo,xxOO)
    fa = zeros(size(ha))
    R = 1:size(xxoo)[1]
    rocca = 1:size(xxoo)[4]
    @tensor begin
        fa[p,q] = ha[p,q] + xxoo[p,q,k,k] - xoxo[p,k,q,k] + xxOO[p,q,k,k]
    end
    #for p in R
    #    for q in R
    #        tsum = 0
    #        for k in rocca
    #            tsum += xxoo[p,q,k,k] - xoxo[p,k,q,k]
    #        end
    #        for k in roccb
    #            tsum += XXOO[p,q,k,k]
    #        end
    #        fa[p,q] += tsum
    #    end
    #end
    return fa
end
function form_fb(hb,XXOO,XOXO,xxoo)
    fb = zeros(size(hb))
    @tensor begin
        fb[p,q] = hb[p,q] + XXOO[p,q,k,k] - XOXO[p,k,q,k] + xxoo[p,q,k,k]
    end
    return fb
end
function form_Fae(fa_vv,fa_ov,oovv,oOvV,vovv,vOvV,tia,tIA,tiatia,tiatIA,tijab,tiJaB)
    Fae = deepcopy(fa_vv) - diagm(diag(fa_vv))
    @tensor begin
        Fae[a,e] += (-fa_ov[m,e]*tia[m,a]
                     +tia[m,f]*(vovv[a,m,e,f] - vovv[a,m,f,e])
                    +tIA[m,f]*vOvV[a,m,e,f]
                    -(1/2)*(tijab[m,n,a,f] + (1/2)*(tiatia[m,n,a,f] - tiatia[m,n,f,a]))*(oovv[m,n,e,f] - oovv[m,n,f,e])
                    -(tiJaB[m,n,a,f] + (1/2)*tiatIA[m,n,a,f])*oOvV[m,n,e,f])
    end
    return Fae
end
function form_FAE(fb_vv,fb_ov,OOVV,OoVv,VOVV,VoVv,tia,tIA,tIAtIA,tIAtia,tIJAB,tiJaB)
    FAE = deepcopy(fb_vv) - diagm(diag(fb_vv))
    @tensor begin
        FAE[a,e] += (-fb_ov[m,e]*tIA[m,a]
                     +tIA[m,f]*(VOVV[a,m,e,f] - VOVV[a,m,f,e])
                     +tia[m,f]*VoVv[a,m,e,f]
                     -(1/2)*(tIJAB[m,n,a,f]+(1/2)*(tIAtIA[m,n,a,f] - tIAtIA[m,n,f,a]))*(OOVV[m,n,e,f] - OOVV[m,n,f,e])
                     -(tiJaB[n,m,f,a]+(1/2)*tIAtia[m,n,a,f])*OoVv[m,n,e,f])
    end
    return FAE
end
function form_Fmi(fa_oo,fa_ov,ooov,oovo,oOoV,oovv,oOvV,tia,tIA,tiatia,tiatIA,tijab,tiJaB)
    Fmi = deepcopy(fa_oo) - diagm(diag(fa_oo))
    @tensor begin
        Fmi[m,i] += ((1/2)*fa_ov[m,e]*tia[i,e]
                     +tia[n,e]*(ooov[m,n,i,e] - oovo[m,n,e,i])
                     +tIA[n,e]*oOoV[m,n,i,e]
                     +(1/2)*(tijab[i,n,e,f] + (1/2)*(tiatia[i,n,e,f] - tiatia[i,n,f,e]))*(oovv[m,n,e,f] - oovv[m,n,f,e])
                     +(tiJaB[i,n,e,f] + (1/2)*tiatIA[i,n,e,f])*oOvV[m,n,e,f])
    end
    return Fmi
end
function form_FMI(fb_oo,fb_ov,OOOV,OOVO,OoOv,OOVV,OoVv,tia,tIA,tIAtIA,tIAtia,tIJAB,tiJaB)
    FMI = deepcopy(fb_oo) - diagm(diag(fb_oo))
    @tensor begin
        FMI[m,i] += ((1/2)*fb_ov[m,e]*tIA[i,e]
                     +tIA[n,e]*(OOOV[m,n,i,e] - OOVO[m,n,e,i])
                     +tia[n,e]*OoOv[m,n,i,e]
                     +(1/2)*(tIJAB[i,n,e,f] + (1/2)*(tIAtIA[i,n,e,f] - tIAtIA[i,n,f,e]))*(OOVV[m,n,e,f] - OOVV[m,n,f,e])
                     +(tiJaB[n,i,f,e] + (1/2)*tIAtia[i,n,e,f])*OoVv[m,n,e,f])
    end
    return FMI
end
function form_Fme(fa,oovv,oOvV,tia,tIA)
end
function form_FME(fb,OOVV,OoVv,tia,tIA)
end
function form_Wmnij(oooo,ooov,oovv,tia,tiatia,tijab)
end
function form_WmNiJ(oOoO,oOoV,oOOv,oOvV,oOVv,tia,tIA,tiatIA,tiJaB,tiJAb)
end
end #module
