"""
    Fermi.CoupledCluster.RCCSD

## Functions

    Fermi.CoupledCluster.RCCSD.do_rccsd
"""
module RCCSD
using Fermi.Wavefunction
using Fermi.Transformation
using Fermi.Output
using Printf
using Fermi
using TensorOperations
include("Denominators.jl")

"""
    Fermi.CoupledCluster.RCCSD.do_rccsd(refWfn::Wfn; kwargs...)

## Arguments
    refWfn::Wfn
reference wavefunction to compute RCCSD on.

    kwargs::Any
keyword arguments for various options.

## Kwargs
    :doprint::Bool  print or no? deprecated
    :maxit::Int     max number of iterations for RCCSD equations
    :return_T::Bool do return the converged cluster amplitudes?
    :diis::Bool     do DIIS extrapolation? currently dummy  
"""
function do_rccsd(refWfn::Wfn; kwargs...)
    Fermi.CoupledCluster.print_header()
    @output "*   executing CCSD\n"
    maxit = 40
    return_T = false
    diis = false
    #defaults = Dict(
    #                :maxit => 40,
    #                :return_T => false,
    #                :diis => false
    #               )
    #for arg in (:maxit,:return_T,:diis)
    #    if arg in keys(kwargs)
    #        @eval $arg = $(kwargs[arg])
    #    else
    #        @eval $arg = $(defaults[arg])
    #    end
    #end
    set_zero_subnormals(true)
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    uvsr = refWfn.ao_eri
    Cav = refWfn.Cav
    Cao = refWfn.Cao
    dtt = eltype(uvsr)

    @output "    Forming MO basis integrals ... "
    t = @elapsed begin
        vvvv, ovvv, vovv,
        vvov, vvvo, oovv,
        ovvo, vovo, ovov,
        voov, ooov, oovo,
        ovoo, vooo, oooo = make_rccsd_integrals(uvsr,Cao,Cav)
    end
    @output "done in {:>5.2f}s\n" t
    epsa = refWfn.epsa
    T2 = zeros(dtt, nocc, nocc, nvir, nvir)
    T1 = zeros(dtt,nocc,nvir)
    Fae = zeros(nvir,nvir)
    Fmi = zeros(nocc,nocc)
    Fme = zeros(nocc,nvir)
    Wabef = zeros(nvir,nvir,nvir,nvir)
    Wmnij = zeros(nocc,nocc,nocc,nocc)
    WmBeJ = zeros(nocc,nvir,nvir,nocc)
    WmBEj = zeros(nocc,nvir,nvir,nocc)
    Dia = form_Dia(T1, epsa)
    Dijab = form_Dijab(T2, epsa)
    T2 .= oovv ./ Dijab
    @tensor begin
        tiatia[m,n,a,f] := T1[m,a]*T1[n,f]
    end
    @output "    @MP2 {:>20.17f}\n" ccenergy(oovv,T1,tiatia,T2)
    @output "    Starting CC iterations\n\n"
    @output "{:>7s} {:>20}\n" "Iter." "E[CCSD]"
    @output repeat("~",80)*"\n"
    for i in 1:maxit
        _T1,_T2 = cciter(Fae,Fmi,Fme,
               Wabef,Wmnij,WmBeJ,WmBEj,
               oovv,     vovv,
               ovvv,ooov,oooo,
               oovo,vvvv,ovvo,
                         voov,
               vovo,ovov,
               vvvo,vvov,ovoo,
               vooo,
               T1,tiatia,T2,Dia,Dijab)
        T1 = _T1
        T2 = _T2
        @tensoropt begin
            tiatia[m,n,a,f] = T1[m,a]*T1[n,f]
        end
        output("    {:<3d} {:>20.17f}\n",i,ccenergy(oovv,T1,tiatia,T2))

    end
    output(repeat("~",80)*"\n")
    @output "    CC iterations completed\n"
    @output "    CC final energy:\n"
    @output "        @CCSD {:>20.17f}\n" ccenergy(oovv,T1,tiatia,T2)
    if return_T
        return ccenergy(oovv,T1,tiatia,T2),T2
    else
        return ccenergy(oovv,T1,tiatia,T2)
    end
end
function make_rccsd_integrals(uvsr,Cao,Cav)
    vvvv = permutedims(tei_transform(uvsr, Cav, Cav, Cav, Cav, "vvvv"),[1,3,2,4])
    ovvv = permutedims(tei_transform(uvsr, Cao, Cav, Cav, Cav, "ovvv"),[1,3,2,4])
    vovv = permutedims(tei_transform(uvsr, Cav, Cav, Cao, Cav, "vovv"),[1,3,2,4])
    vvov = permutedims(tei_transform(uvsr, Cav, Cao, Cav, Cav, "vvov"),[1,3,2,4])
    vvvo = permutedims(tei_transform(uvsr, Cav, Cav, Cav, Cao, "vvvo"),[1,3,2,4])
    oovv = permutedims(tei_transform(uvsr, Cao, Cav, Cao, Cav, "oovv"),[1,3,2,4])
    ovvo = permutedims(tei_transform(uvsr, Cao, Cav, Cav, Cao, "ovvo"),[1,3,2,4])
    vovo = permutedims(tei_transform(uvsr, Cav, Cav, Cao, Cao, "vovo"),[1,3,2,4])
    ovov = permutedims(tei_transform(uvsr, Cao, Cao, Cav, Cav, "ovov"),[1,3,2,4])
    voov = permutedims(tei_transform(uvsr, Cav, Cao, Cao, Cav, "voov"),[1,3,2,4])
    ooov = permutedims(tei_transform(uvsr, Cao, Cao, Cao, Cav, "ooov"),[1,3,2,4]) 
    oovo = permutedims(tei_transform(uvsr, Cao, Cav, Cao, Cao, "oovo"),[1,3,2,4])
    ovoo = permutedims(tei_transform(uvsr, Cao, Cao, Cav, Cao, "ovoo"),[1,3,2,4])
    vooo = permutedims(tei_transform(uvsr, Cav, Cao, Cao, Cao, "vooo"),[1,3,2,4])
    oooo = permutedims(tei_transform(uvsr, Cao, Cao, Cao, Cao, "oooo"),[1,3,2,4])
    vvov = permutedims(vvov,[4,1,2,3])
    vvvo = permutedims(vvvo,[3,1,2,4])
    vovv = permutedims(vovv,[2,1,3,4])
    vooo = permutedims(vooo,[2,1,3,4])
    return vvvv, ovvv, vovv,
           vvov, vvvo, oovv,
           ovvo, vovo, ovov,
           voov, ooov, oovo,
           ovoo, vooo, oooo
end
function ccenergy(ijab,tia,tiatia,tijab)
    @tensoropt begin
        ecc[] := ijab[i,j,a,b]*(2*tijab[i,j,a,b] + 2*tiatia[i,j,a,b] 
                              - tijab[j,i,a,b] - tiatia[j,i,a,b])
    end
    return ecc[]
end
function cciter(Fae,Fmi,Fme,
                Wabef,Wmnij,WmBeJ,WmBEj,
                oovv,     vovv,
                ovvv,ooov,oooo,
                oovo,vvvv,ovvo,
                          voov,
                vovo,ovov,
                vvvo,vvov,ovoo,
                vooo,
                tia,tiatia,tijab,Dia,Dijab)
    form_Fae!(Fae,ovvv,vovv,oovv,tia,tiatia,tijab)
    form_Fmi!(Fmi,ooov,oovv,tia,tiatia,tijab)
    form_Fme!(Fme,oovv,tia)
    form_Wabef!(Wabef,vvvv,vovv,ovvv,oovv,tia,tiatia,tijab)
    form_Wmnij!(Wmnij,oooo,ooov,oovo,oovv,tia,tiatia,tijab)
    form_WmBeJ!(WmBeJ,ovvo,ovvv,oovo,oovv,tia,tiatia,tijab)
    form_WmBEj!(WmBEj,vovo,ovvv,oovo,oovv,tia,tiatia,tijab)
    _tia = form_T1(Fae,Fmi,Fme,voov,ovov,ooov,vovv,tia,tijab,Dia)
    _tijab = form_T2(Fae,Fmi,Fme,Wabef,Wmnij,WmBeJ,WmBEj,oovv,
                     ovvo,vovo,vvvo,vvov,ovoo,vooo,tia,tiatia,tijab,Dijab)
    #tia = nothing
    #tijab = nothing
    return _tia,_tijab
end
                
function form_Dia(T1, epsa)
    Dia = zeros(size(T1))
    nocc = size(T1)[1]
    nvir = size(T1)[2]
    for i in 1:nocc
        for a in 1:nvir
            aa = a + nocc
            Dia[i,a] = epsa[i] - epsa[aa]
        end
    end
    return Dia
end
function form_Fae!(Fae,maef,amef,mnef,tia,tiatia,tijab)
    #Fae .= 0.0
    @tensoropt begin
        Fae[a,e] = (tia[m,f]*(2*amef[m,a,e,f] - maef[m,a,e,f])
                    - ((tijab[m,n,a,f] + 0.5*tiatia[m,n,a,f])*(2*mnef[m,n,e,f] - mnef[n,m,e,f])))
    end
    return Fae
end
function form_Fmi!(Fmi,mnie,mnef,tia,tiatia,tijab)
    #Fmi .= 0.0
    @tensoropt begin
        Fmi[m,i] = (tia[n,e]*(2*mnie[m,n,i,e] - mnie[n,m,i,e])
                    + ((tijab[i,n,e,f] + 0.5*tiatia[i,n,e,f])
                       *(2*mnef[m,n,e,f] - mnef[m,n,f,e])))
    end
    return Fmi
end
function form_Fme!(Fme,mnef,tia)
    #Fme .= 0.0
    @tensoropt begin
        Fme[m,e] = tia[n,f]*(2*mnef[m,n,e,f] - mnef[n,m,e,f])
    end
    return Fme
end
function form_Wmnij!(Wmnij,mnij,mnie,mnej,mnef,tia,tiatia,tijab)
    #Wmnij .= 0.0
    @tensoropt begin
        Wmnij[m,n,i,j] = (mnij[m,n,i,j] + tia[j,e]*mnie[m,n,i,e]
                          + tia[i,e]*mnej[m,n,e,j]
                          + 0.5*(tijab[i,j,e,f] + tiatia[i,j,e,f])*mnef[m,n,e,f])
    end
    return Wmnij
end
function form_Wabef!(Wabef,abef,amef,mbef,mnef,tia,tiatia,tijab)
    #Wabef .= 0.0
    @tensoropt begin
        Wabef[a,b,e,f] = (abef[a,b,e,f] - tia[m,b]*amef[m,a,e,f]
                          - tia[m,a]*mbef[m,b,e,f]
                          + 0.5*(tijab[m,n,a,b] + tiatia[m,n,a,b])*mnef[m,n,e,f])
    end
    return Wabef
end
function form_WmBeJ!(WmBeJ,mbej,mbef,mnej,mnef,tia,tiatia,tijab)
    #WmBeJ .= 0.0
    @tensoropt begin
        WmBeJ[m,b,e,j] = (mbej[m,b,e,j] + tia[j,f]*mbef[m,b,e,f]
                          - tia[n,b]*mnej[m,n,e,j]
                          - 0.5*(tijab[j,n,f,b] + 2*tiatia[j,n,f,b])*mnef[m,n,e,f]
                          + 0.5*tijab[n,j,f,b]*(2*mnef[m,n,e,f] - mnef[n,m,e,f]))
    end
    return WmBeJ
end
function form_WmBEj!(WmBEj,bmej,mbfe,nmej,nmef,tia,tiatia,tijab)
    #WmBEj .= 0.0
    @tensoropt begin
        WmBEj[m,b,e,j] = (-1*bmej[b,m,e,j] - tia[j,f]*mbfe[m,b,f,e]
                          + tia[n,b]*nmej[n,m,e,j]
                          + (0.5*tijab[j,n,f,b] + tiatia[j,n,f,b])*nmef[n,m,e,f])
    end
    return WmBEj
end
function form_T1(Fae,Fmi,Fme,amie,maie,mnie,amef,tia,tijab,Dia)
    #_tia = zeros(size(tia))
    @tensoropt begin
        _tia[i,a] := (tia[i,e]*Fae[a,e] - tia[m,a]*Fmi[m,i]
                     + Fme[m,e]*(2*tijab[i,m,a,e] - tijab[m,i,a,e])
                     + tia[m,e]*(2*amie[a,m,i,e] - maie[m,a,i,e])
                     - tijab[m,n,a,e]*(2*mnie[m,n,i,e] - mnie[n,m,i,e])
                     + tijab[i,m,e,f]*(2*amef[m,a,e,f] - amef[m,a,f,e]))
    end
    _tia .= _tia ./ Dia
    return _tia
end
function form_T2(Fae,Fmi,Fme,Wabef,Wmnij,WmBeJ,WmBEj,ijab,mbej,amej,abej,abie,mbij,amij,tia,tiatia,tijab,Dijab)
    #_tijab = zeros(size(tijab))
    _tia = permutedims(tia,[2,1])
    @tensoropt (i=>x,j=>x,m=>x,n=>x,a=>11*x,b=>11*x,e=>11*x,f=>11*x) begin
        _tijab[i,j,a,b] := (ijab[i,j,a,b] 
                           + tijab[i,j,a,e]*(Fae[b,e] - 0.5*tia[m,b]*Fme[m,e])
                           + tijab[i,j,e,b]*(Fae[a,e] - 0.5*tia[m,a]*Fme[m,e])
                           - tijab[i,m,a,b]*(Fmi[m,j] + 0.5*tia[j,e]*Fme[m,e])
                           - tijab[m,j,a,b]*(Fmi[m,i] + 0.5*tia[i,e]*Fme[m,e])
                           + (tijab[m,n,a,b] + tiatia[m,n,a,b])*Wmnij[m,n,i,j]
                           + (tijab[i,j,e,f] + tiatia[i,j,e,f])*Wabef[a,b,e,f]
                           + ((tijab[i,m,a,e] - tijab[m,i,a,e])*WmBeJ[m,b,e,j]
                              - tiatia[i,m,e,a]*mbej[m,b,e,j])
                           + tijab[i,m,a,e]*(WmBeJ[m,b,e,j] + WmBEj[m,b,e,j])
                           + (tijab[m,i,b,e]*WmBEj[m,a,e,j] 
                              - tiatia[i,m,e,b]*amej[a,m,e,j])
                           + (tijab[m,j,a,e]*WmBEj[m,b,e,i]
                              - tiatia[j,m,e,a]*amej[b,m,e,i])
                           + ((tijab[j,m,b,e] - tijab[m,j,b,e])*WmBeJ[m,a,e,i]
                              - tiatia[j,m,e,b]*mbej[m,a,e,i])
                           + tijab[j,m,b,e]*(WmBeJ[m,a,e,i] + WmBEj[m,a,e,i])
                           + _tia[e,i]*abej[e,a,b,j]
                           + _tia[e,j]*abie[e,a,b,i]
                           - tia[m,a]*mbij[m,b,i,j]
                           - tia[m,b]*amij[m,a,i,j]
                          )
    end
    _tijab .= _tijab ./ Dijab
    return _tijab
end
end #module
