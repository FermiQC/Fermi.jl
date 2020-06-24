module RCCD2
using LoopVectorization
using Base.Threads
using SparseArrays
using Fermi.Wavefunction
using TensorOperations
using Fermi.Transformation
include("Denominators.jl")
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
    T2 = zeros(T, nocc, nocc, nvir, nvir)
    Dijab = form_Dijab(T2, epsa)
    T2_init!(T2, ovov, Dijab)
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
        T2 = cciter(
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
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for i in rocc
        for j in rocc
            #cache = ijab[i,j,:,:]
            for a in rvir
                for b in rvir
                    ecc += ijab[i,j,a,b] * 2 * tiJaB[i, j, a, b]
                    ecc -= ijab[i,j,a,b] * tiJaB[j, i, a, b]
                end
            end
        end
    end
    return ecc
end

@inbounds @fastmath function cciter(
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

    tiJaB_d = form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, oovv, Dijab)
    return tiJaB_d
end

function T2_init!(tiJaB, iajb, Dijab)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    for i in rocc
        for j in rocc
            cache = iajb[i,:,j,:]
            for a in rvir
                for b in rvir
                    tiJaB[i,j,a,b] = cache[a,b] / Dijab[i,j,a,b]
                end
            end
        end
    end
end

function form_Fae(tiJaB, menf)
    dt = eltype(menf)
    nvir = size(tiJaB, 4)
    Fae = zeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, menf)
    return Fae
end
function form_Fae!(Fae, tiJaB, mnef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    Fae .= 0.0
    #cache1 = zeros(eltype(menf), nocc, nocc)
    #cache2 = zeros(eltype(menf), nocc, nocc)
    #mnef = permutedims(menf,[1,3,2,4])
    @tensoropt begin
        Fae[a,e] = -1*(mnef[m,n,e,f]*(2*tiJaB[m,n,a,f] - tiJaB[n,m,a,f]))
    end
    #for f in rvir
    #    for a in rvir
    #        _tiJaB = tiJaB[:,:,a,f] .* 2.0
    #        _tiJaBt = permutedims(tiJaB[:,:,a,f],[2,1])
    #        for e in rvir
    #            cache = menf[:,e,:,f]
    #            @views Fae[a, e] -= dot(cache[:,:], _tiJaB[:,:] .- _tiJaBt[:,:])
    #            #for n in rocc
    #            #    Fae[a, e] -= 
    #            #    @simd for m in rocc
    #            #        Fae[a, e] -=
    #            #            cache[m,n] * (_tiJaB[m,n] - _tiJaBt[m,n])
    #            #    end
    #            #end
    #        end
    #    end
    #end
end

function form_Fmi(tiJaB, menf)
    dt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Fmi = zeros(dt, nocc, nocc)
    form_Fmi!(Fmi, tiJaB, menf)
    return Fmi
end
function form_Fmi!(Fmi, tiJaB, mnef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    Fmi .= 0.0
    #mnef = permutedims(menf,[1,3,2,4])
    @tensoropt begin
        Fmi[m,i] = mnef[m,n,e,f]*(2*tiJaB[i,n,e,f] - tiJaB[i,n,f,e])
    end
    #for f in rvir
    #    for e in rvir
    #        for i in rocc
    #            cache = menf[:,e,:,f]
    #            for n in rocc
    #                for m in rocc
    #                    Fmi[m, i] +=
    #                        cache[m,n] * (2 * tiJaB[i, n, e, f] - tiJaB[i, n, f, e])
    #                end
    #            end
    #        end
    #    end
    #end
    return Fmi
end

function form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, ijab, Dijab)
    dtt = eltype(tiJaB_i)
    nocc = size(Wmnij, 1)
    nvir = size(tiJaB_i, 4)
    rocc = 1:nocc
    rvir = 1:nvir 
    #ijab = permutedims(iajb,[1,3,2,4])
    tiJaB_d = zeros(size(ijab[:,:,:,:]))
    #tiJaB_d .= ijab
    for b in rvir
        for a in 1:nvir
            for j in rocc
                for i in j:nocc
                    tiJaB_d[i,j,a,b] = ijab[i,j,a,b]
                    tiJaB_d[j,i,a,b] = ijab[j,i,a,b]
                end
            end
        end
    end

    tiJaB_d_temp = zeros(size(tiJaB_d))
    tiJaB_d_temp2 = zeros(size(tiJaB_d))
    tiJaB_d_temp3 = zeros(size(tiJaB_d))
    tiJaB_d_temp4 = zeros(size(tiJaB_d))
    tiJaB_d_temp5 = zeros(size(tiJaB_d))
    tiJaB_d_temp6 = zeros(size(tiJaB_d))
    tiJaB_i_2 = permutedims(tiJaB_i,(4,1,2,3))
    WmBeJ_2 = permutedims(WmBeJ,(3,1,2,4))
    WmBEj_2 = permutedims(WmBEj,(3,1,2,4))
    #Wabef = permutedims(Wabef,(4,3,2,1))
    tiJaB_2 = tensorcopy(tiJaB_i,(1,2,3,4),(4,3,1,2))
    @inbounds @fastmath Threads.@threads for b in rvir
        @avx for a in b:nvir
            for j in rocc
                for i in rocc
                    temp = 0.0
                    for m in rocc
                        for e in rvir
                            temp += tiJaB_i_2[e,i,m,a]*WmBeJ_2[e,m,b,j]*2
                            temp -= tiJaB_i_2[e,m,i,a]*WmBeJ_2[e,m,b,j]
                            temp += tiJaB_i_2[e,i,m,a]*WmBEj_2[e,m,b,j]
                            temp += tiJaB_i_2[e,m,i,b]*WmBEj_2[e,m,a,j]
                            temp += tiJaB_i_2[e,m,j,a]*WmBEj_2[e,m,b,i]
                            temp += tiJaB_i_2[e,j,m,b]*WmBeJ_2[e,m,a,i]*2
                            temp -= tiJaB_i_2[e,m,j,b]*WmBeJ_2[e,m,a,i]
                            temp += tiJaB_i_2[e,j,m,b]*WmBEj_2[e,m,a,i]
                        end
                    end
                    tiJaB_d_temp[i,j,a,b] = temp
                    tiJaB_d_temp[j,i,b,a] = temp
                end
            end
        end
    end
    println(1)
    #@inbounds @fastmath Threads.@threads for b in rvir
    #    @avx for a in rvir
    #        for j in rocc
    #            for i in rocc
    #                temp = 0.0
    #                for m in rocc
    #                    for e in rvir
    #                        temp += tiJaB_i_2[e,m,j,a]*WmBEj_2[e,m,b,i]
    #                        temp += tiJaB_i_2[e,j,m,b]*WmBeJ_2[e,m,a,i]*2
    #                        temp -= tiJaB_i_2[e,m,j,b]*WmBeJ_2[e,m,a,i]
    #                        temp += tiJaB_i_2[e,j,m,b]*WmBEj_2[e,m,a,i]
    #                    end
    #                end
    #                tiJaB_d_temp2[i,j,a,b] = temp
    #            end
    #        end
    #    end
    #end
    @inbounds @fastmath Threads.@threads for b in rvir
        @avx for a in b:nvir
            for j in rocc
                for i in rocc
                    temp = zero(eltype(tiJaB_i))
                    for n in rocc
                        for m in rocc
                            temp += tiJaB_i[m,n,a,b]*Wmnij[m,n,i,j]
                        end
                    end
                    tiJaB_d_temp3[i,j,a,b] = temp
                    tiJaB_d_temp3[j,i,b,a] = temp
                end
            end
        end
    end
    println(2)
    @inbounds @fastmath Threads.@threads for b in rvir
        @avx for a in b:nvir
            for j in rocc
                for i in rocc
                    temp = zero(eltype(tiJaB_i))
                    for m in rocc
                        temp -= tiJaB_i[i,m,a,b]*Fmi[m,j]
                        temp -= tiJaB_i[m,j,a,b]*Fmi[m,i]
                    end
                    tiJaB_d_temp4[i,j,a,b] = temp
                    tiJaB_d_temp4[j,i,b,a] = temp
                end
            end
        end
    end
    println(3)
    #tiJaB_2 = permutedims(tiJaB_i,(4,3,1,2))
    #@tensor tiJaB_d_temp5[i,j,a,b] = tiJaB_2[f,e,i,j]*Wabef[a,b,e,f]
    @inbounds @fastmath Threads.@threads for b in rvir
        @avx for a in b:nvir
            for j in rocc
                for i in rocc
                    temp = zero(eltype(tiJaB_i))
                    for e in rvir
                        for f in rvir
                            temp += tiJaB_2[f,e,i,j]*Wabef[a,b,e,f]
                        end
                    end
                    tiJaB_d_temp5[i,j,a,b] = temp
                    tiJaB_d_temp5[j,i,b,a] = temp
                end
            end
        end
    end
    println(4)
    #@tensor tiJaB_d_temp6[i,j,a,b] = tiJaB_i[i,j,a,e]*Fae[b,e] + tiJaB_i[j,i,b,e]*Fae[a,e]
    @inbounds @fastmath Threads.@threads for b in rvir
        @avx for a in b:nvir
            for j in rocc
                for i in rocc
                    temp = zero(eltype(tiJaB_i))
                    for e in rvir
                        temp += tiJaB_i[i,j,a,e]*Fae[b,e]
                        temp += tiJaB_i[j,i,b,e]*Fae[a,e]
                    end
                    tiJaB_d_temp6[i,j,a,b] = temp
                    tiJaB_d_temp6[j,i,b,a] = temp
                end
            end
        end
    end
    println(5)
    tiJaB_d .+= tiJaB_d_temp
    tiJaB_d .+= tiJaB_d_temp2
    tiJaB_d .+= tiJaB_d_temp3
    tiJaB_d .+= tiJaB_d_temp4
    tiJaB_d .+= tiJaB_d_temp5
    tiJaB_d .+= tiJaB_d_temp6
    #@tensoropt begin
    #    tiJaB_d[i,j,a,b] = (ijab[i,j,a,b] + tiJaB_i[i,j,a,e]*Fae[b,e] + tiJaB_i[j,i,b,e]*Fae[a,e]
    #                       - tiJaB_i[i,m,a,b]*Fmi[m,j] - tiJaB_i[m,j,a,b]*Fmi[m,i]
    #                       + tiJaB_i[m,n,a,b]*Wmnij[m,n,i,j]
    #                       + tiJaB_i[i,j,e,f]*Wabef[a,b,e,f] 
    #                       + tiJaB_i[i,m,a,e]*WmBeJ[m,b,e,j]*2
    #                       - tiJaB_i[m,i,a,e]*WmBeJ[m,b,e,j]
    #                       + tiJaB_i[i,m,a,e]*WmBEj[m,b,e,j]
    #                       + tiJaB_i[m,i,b,e]*WmBEj[m,a,e,j]
    #                       + tiJaB_i[m,j,a,e]*WmBEj[m,b,e,i]
    #                       + tiJaB_i[j,m,b,e]*WmBeJ[m,a,e,i]*2
    #                       - tiJaB_i[m,j,b,e]*WmBeJ[m,a,e,i]
    #                       + tiJaB_i[j,m,b,e]*WmBEj[m,a,e,i])

    #end
    #for b in rvir
    #    for a in rvir
    #        cache_iajb = iajb[:,a,:,b]
    #        cache_tijab = tiJaB_d[:,:,a,b]
    #        for j in rocc
    #            for i in rocc
    #                temp = 0.0
    #                temp += ijab[i,j,a,b]#cache_iajb[i, j]
    #                #@views temp += dot(tiJaB_i[i,j,:,:], Wabef[a,b,:,:])
    #                @views temp += dot(tiJaB_i[i,j,a,:], Fae[b,:])
    #                @views temp += dot(tiJaB_i[j,i,b,:], Fae[a,:])
    #                #for e in rvir
    #                #    temp += tiJaB_i[i, j, a, e] * Fae[b, e]
    #                #    temp += tiJaB_i[j, i, b, e] * Fae[a, e]
    #                #    for f in rvir
    #                #        temp += tiJaB_i[i, j, e, f] * Wabef[a, b, e, f]
    #                #    end
    #                #end
    #                for m in rocc
    #                    temp -= tiJaB_i[i, m, a, b] * Fmi[m, j]
    #                    temp -= tiJaB_i[m, j, a, b] * Fmi[m, i]
    #                    for n in rocc
    #                        temp += tiJaB_i[m, n, a, b] * Wmnij[m, n, i, j]
    #                    end
    #                    for e in rvir
    #                        c1 = WmBeJ[m,b,e,j]
    #                        c2 = tiJaB_i[i,m,a,e]
    #                        c3 = WmBeJ[m,a,e,i]
    #                        c4 = tiJaB_i[j,m,b,e]
    #                        temp += (2 *( c2 * c1 + ( c4 * c3)) - (tiJaB_i[m, i, a, e] * c1) + (c2 * WmBEj[m, b, e, j])
    #                                 + (tiJaB_i[m, i, b, e] * WmBEj[m, a, e, j])
    #                                 + (tiJaB_i[m, j, a, e] * WmBEj[m, b, e, i])
    #                                 - (tiJaB_i[m, j, b, e] * c3)
    #                                 + (c4 * WmBEj[m, a, e, i]))
    #                    end
    #                end
    #                cache_tijab[i,j] += temp
    #            end
    #        end
    #        tiJaB_d[:,:,a,b] += cache_tijab# ./ Dijab[:,:,a,b]
    #    end
    #end
    #tiJaB_i = nothing
    #tiJaB_d = convert(Array,tiJaB_d)
    tiJaB_d .= tiJaB_d ./ Dijab
    return tiJaB_d
end
function form_Wmnij(minj, menf, tiJaB)
    dtt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Wmnij = zeros(dtt, nocc, nocc, nocc, nocc)
    form_Wmnij!(Wmnij, minj, menf, tiJaB)
    return Wmnij
end
function form_Wmnij!(Wmnij, mnij, mnef, tiJaB)
    eps = 1E-10
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    #Wmnij = convert(SharedArray,Wmnij)
    #mnij = permutedims(minj,[1,3,2,4])
    #mnef = permutedims(menf,[1,3,2,4])
    #tiJaB2 = zeros(nocc^2,nvir^2)
    #mnef2 = zeros(nocc^2,nvir^2)
    #for i=1:nocc,j=1:nocc,a=1:nvir,b=1:nvir
    #    tiJaB2[(i-1)*nocc+j,(a-1)*nvir+b] = abs(tiJaB[i,j,a,b]) > eps ? tiJaB[i,j,a,b] : 0
    #    mnef2[(i-1)*nocc+j,(a-1)*nvir+b]  = abs(mnef[i,j,a,b]) > eps ? mnef[i,j,a,b] : 0
    #end
    #print(tiJaB2)
    #test = sparse(tiJaB2)
    #test2 = sparse(mnef2)
    #dropzeros!(test)
    #dropzeros!(test2)
    #print(size(test))
    #print(nnz(test))
    #print(size(test2))
    #print(nnz(test2))
    #tiJaBmnef2 = test*transpose(mnef2)
    #tiJaBmnef  = zeros(nocc,nocc,nocc,nocc)
    #for m=1:nocc,n=1:nocc,i=1:nocc,j=1:nocc
    #    tiJaBmnef[m,n,i,j] = tiJaBmnef2[(i-1)*nocc+j,(m-1)*nocc+n]
    #end
    @tensoropt begin
        Wmnij[m,n,i,j] = mnij[m,n,i,j] + tiJaB[i,j,e,f]*mnef[m,n,e,f]/2
    end
    #for m in rocc
    #    for n in rocc
    #        cache_minj = minj[m,:,n,:]
    #        cache_menf = menf[m,:,n,:]
    #        for i in rocc
    #            for j in rocc
    #                Wmnij[m, n, i, j] = cache_minj[i, j]
    #                @views Wmnij[m, n, i, j] += dot(tiJaB[i,j,:,:],cache_menf[:,:])/2.0
    #                #for f in rvir
    #                #    @simd for e in rvir
    #                #        Wmnij[m, n, i, j] += tiJaB[i, j, e, f] * cache_menf[e, f]/2.0 
    #                #    end
    #                #end
    #            end
    #        end
    #    end
    #end
    #Wmnij = convert(Array,Wmnij)
end

function form_Wabef(aebf, mnef, tiJaB)
    dt = eltype(aebf)
    nvir = size(tiJaB, 4)
    Wabef = zeros(dt, nvir, nvir, nvir, nvir)
    form_Wabef!(Wabef, aebf, mnef, tiJaB)
    return Wabef
end
function form_Wabef!(Wabef, abef, mnef, tiJaB)
    dtt = eltype(Wabef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    #mnef = permutedims(menf,[1,3,2,4])
    #abef = permutedims(aebf,[1,3,2,4])
    @tensoropt begin
        Wabef[a,b,e,f] = tiJaB[m,n,a,b]*mnef[m,n,e,f]/2 + abef[a,b,e,f]
    end
    #Threads.@threads for f in rvir
    #    @avx for e in rvir
    #        for b in rvir
    #            for a in rvir
    #                _Wabef = zero(eltype(Wabef))
    #                for n in rocc
    #                    for m in rocc
    #                        _Wabef += tiJaB[m,n,a,b]*mnef[m,n,e,f] #+ abef[a,b,e,f]
    #                    end
    #                end
    #                Wabef[a,b,e,f] = _Wabef/2
    #            end
    #        end
    #    end
    #end
    #Wabef .+= abef
    #for f in rvir
    #    for e in rvir
    #        cache_aebf = aebf[:,e,:,f]
    #        cache_menf = menf[:,e,:,f]
    #        for b in rvir
    #            for a in rvir
    #                #Wabef[a, b, e, f] = 0 
    #                @views Wabef[a, b, e, f] = dot(tiJaB[:,:,a,b], cache_menf[:,:])/2.0
    #                #for n in rocc
    #                #    @simd for m in rocc
    #                #        Wabef[a, b, e, f] += tiJaB[m, n, a, b] * cache_menf[m, n] 
    #                #    end
    #                #end
    #                #Wabef[a, b, e, f] /= 2.0
    #                Wabef[a, b, e, f] += cache_aebf[a, b]
    #            end
    #        end
    #    end
    #end
end

function form_WmBeJ(mebj, iajb, tiJaB)
    dtt = eltype(iajb)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBeJ = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBeJ!(WmBeJ, mebj, iajb, tiJaB)
    return WmBeJ
end
function form_WmBeJ!(WmBeJ, mbej, mnef, tiJaB)
    dtt = eltype(WmBeJ)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    #WmBeJ .= 0 
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    #mbej = permutedims(mebj,[1,3,2,4])
    #mnef = permutedims(iajb,[1,3,2,4])
    @tensoropt begin
        WmBeJ[m,b,e,j] = (mbej[m,b,e,j] + mnef[m,n,e,f]*(2*tiJaB[n,j,f,b] - tiJaB[j,n,f,b])/2
                          - mnef[n,m,e,f]*tiJaB[n,j,f,b]/2)
    end
    #for e in rvir
    #    for m in rocc
    #        cache_mebj = mebj[m,e,:,:]
    #        cache_iajb1 = iajb[m,e,:,:]
    #        cache_iajb2 = iajb[:,e,m,:]
    #        for b in rvir
    #            for j in rocc
    #                WmBeJ[m, b, e, j] = 0#cache_mebj[b, j]
    #                for f in rvir
    #                    @simd for n in rocc
    #                        WmBeJ[m, b, e, j] +=
    #                            cache_iajb1[n, f] *
    #                            (2 * tiJaB[n, j, f, b] - tiJaB[j, n, f, b])
    #                        WmBeJ[m, b, e, j] -= cache_iajb2[n, f] * tiJaB[n, j, f, b]
    #                    end
    #                end
    #                WmBeJ[m, b, e, j] /= 2.0
    #                WmBeJ[m, b, e, j] += cache_mebj[b, j]
    #            end
    #        end
    #    end
    #end
end


function form_WmBEj(nemf, mjbe, tiJaB)
    dtt = eltype(nemf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBEj = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBEj!(WmBEj, nemf, mjbe, tiJaB)
    return WmBEj
end
function form_WmBEj!(WmBEj, nmef, mbje, tiJaB)
    dtt = eltype(WmBEj)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    #mbje = permutedims(mjbe,[1,3,2,4])
    #nmef = permutedims(nemf,[1,3,2,4])
    @tensoropt begin
        WmBEj[m,b,e,j] = -mbje[m,b,j,e] + tiJaB[j,n,f,b]*nmef[n,m,e,f]/2.0
    end
    #for e in rvir
    #    for m in rocc
    #        cache_mjbe = mjbe[m,:,:,e]
    #        cache_nemf = nemf[:,e,m,:]
    #        for b in rvir
    #            for j in rocc
    #                #WmBEj[m, b, e, j] = 0#-cache_mjbe[j, b]
    #                @views WmBEj[m, b, e, j] = dot(tiJaB[j, :, :, b], cache_nemf[:, :])/2.0
    #                #for f in rvir
    #                #    @simd for n in rocc
    #                #        WmBEj[m, b, e, j] += tiJaB[j, n, f, b] * cache_nemf[n, f] 
    #                #    end
    #                #end
    #                #WmBEj[m, b, e, j] /= 2.0
    #                WmBEj[m, b, e, j] += -cache_mjbe[j, b]
    #            end
    #        end
    #    end
    #end
    return WmBEj
end
end #module
