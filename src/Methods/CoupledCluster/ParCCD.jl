using Distributed
using SharedArrays
function do_rccd(refWfn::Wfn)
    #implicit maxit = 40
    return do_rccd(refWfn, 40)
end
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
function do_rccd(refWfn::Wfn, maxit; doprint::Bool=false, return_T2::Bool=false)
    #goes through appropriate steps to do RCCD
    dtm = @elapsed begin
    set_zero_subnormals(true)
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    ovov = refWfn.ijab
    dtt = Float64#eltype(ovov)
    vvvv = @spawnat :any tei_transform(refWfn.uvsr, refWfn.Cav, refWfn.Cav, refWfn.Cav, refWfn.Cav, "vvvv")
    ovvo = @spawnat :any tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cav, refWfn.Cav, refWfn.Cao, "ovvo")
    oovv = @spawnat :any tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cav, refWfn.Cav, "oovv")
    ooov = @spawnat :any tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cao, refWfn.Cav, "ooov")
    oooo = @spawnat :any tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cao, refWfn.Cao, "oooo")
    oooo = fetch(oooo)
    ooov = fetch(ooov)
    oovv = fetch(oovv)
    ovvo = fetch(ovvo)
    vvvv = fetch(vvvv)

    epsa = refWfn.epsa
    T2 = zeros(dtt, nocc, nocc, nvir, nvir)
    Dijab = form_Dijab(T2, epsa)
    T2_init!(T2, ovov, Dijab)
    if doprint
        println("@RMP2 ", ccenergy(T2, ovov))
    end
    Fae = form_Fae(T2, ovov)
    Fmi = form_Fmi(T2, ovov)
    Wmnij = form_Wmnij(oooo, ovov, T2)
    Wabef = form_Wabef(vvvv, ovov, T2)
    WmBeJ = form_WmBeJ(ovvo, ovov, T2)
    WmBEj = form_WmBEj(ovov, oovv, T2)
    dt = @elapsed for i in 0:maxit-1 #TODO: implement RMS check
        T2 = cciter(
            T2,
            ovov,
            vvvv,
            oooo,
            oovv,
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
            println("$i @CCD ", ccenergy(T2, ovov))
        end
        t1 = Dates.Time(Dates.now())
    end
    if doprint
        println("CCD iterations computed in $dt s")
    end
    end
    if doprint println("CCD complete in $dtm s") end
    if return_T2
        return ccenergy(T2, ovov), T2
    else
        return ccenergy(T2, ovov)
    end
end

function ccenergy(tiJaB, iajb)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for i in rocc
        for j in rocc
            cache = iajb[i,:,j,:]
            for a in rvir
                for b in rvir
                    ecc += cache[a,b] * 2 * tiJaB[i, j, a, b]
                    ecc -= cache[a,b] * tiJaB[j, i, a, b]
                end
            end
        end
    end
    return ecc
end

function cciter(
    tiJaB_i,
    ovov,
    vvvv,
    oooo,
    oovv,
    ovvo,
    Dijab,
    Fae,
    Fmi,
    Wabef,
    Wmnij,
    WmBeJ,
    WmBEj,
)
    f = @spawnat :any form_Wabef!(Wabef, deepcopy(vvvv), deepcopy(ovov), deepcopy(tiJaB_i))
    d = @spawnat :any form_WmBeJ!(WmBeJ, deepcopy(ovvo), deepcopy(ovov), deepcopy(tiJaB_i))
    e = @spawnat :any form_WmBEj!(WmBEj, deepcopy(ovov), deepcopy(oovv), deepcopy(tiJaB_i))
    c = @spawnat :any form_Wmnij!(Wmnij, deepcopy(oooo), deepcopy(ovov), deepcopy(tiJaB_i))
    a = @spawnat :any form_Fae!(Fae, deepcopy(tiJaB_i), deepcopy(ovov))
    b = @spawnat :any form_Fmi!(Fmi, deepcopy(tiJaB_i), deepcopy(ovov))
    fetch(a)
    fetch(b)
    fetch(c)
    fetch(d)
    fetch(e)
    fetch(f)
    tiJaB_d = form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, ovov, Dijab)
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
function form_Fae!(Fae, tiJaB, menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    Fae .= 0.0
    cache1 = zeros(eltype(menf), nocc, nocc)
    cache2 = zeros(eltype(menf), nocc, nocc)
    for f in rvir
        for a in rvir
            _tiJaB = tiJaB[:,:,a,f] .* 2.0
            _tiJaBt = permutedims(tiJaB[:,:,a,f],[2,1])
            for e in rvir
                cache = menf[:,e,:,f]
                @views Fae[a, e] -= dot(cache[:,:], _tiJaB[:,:] .- _tiJaBt[:,:])
                #for n in rocc
                #    Fae[a, e] -= 
                #    @simd for m in rocc
                #        Fae[a, e] -=
                #            cache[m,n] * (_tiJaB[m,n] - _tiJaBt[m,n])
                #    end
                #end
            end
        end
    end
end

function form_Fmi(tiJaB, menf)
    dt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Fmi = zeros(dt, nocc, nocc)
    form_Fmi!(Fmi, tiJaB, menf)
    return Fmi
end
function form_Fmi!(Fmi, tiJaB, menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    Fmi .= 0.0
    for f in rvir
        for e in rvir
            for i in rocc
                cache = menf[:,e,:,f]
                for n in rocc
                    for m in rocc
                        Fmi[m, i] +=
                            cache[m,n] * (2 * tiJaB[i, n, e, f] - tiJaB[i, n, f, e])
                    end
                end
            end
        end
    end
    return Fmi
end

function form_Dijab(tiJaB, F)
    dt = eltype(tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    Dijab = zeros(dt, nocc, nocc, nvir, nvir)
    for i in rocc
        for j in rocc
            for a in rvir
                for b in rvir
                    aa = a + nocc
                    bb = b + nocc
                    Dijab[i, j, a, b] = F[i] + F[j] - F[aa] - F[bb]
                end
            end
        end
    end
    return Dijab
end

function form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, iajb, Dijab)
    dtt = eltype(tiJaB_i)
    nocc = size(Wmnij, 1)
    nvir = size(tiJaB_i, 4)
    tiJaB_d = convert(SharedArray,zeros(dtt, nocc, nocc, nvir, nvir))
    rocc = 1:nocc
    rvir = 1:nvir 
    dt = @elapsed begin
    @sync @distributed for b in rvir
        for a in rvir
            cache_iajb = iajb[:,a,:,b]
            cache_tijab = tiJaB_d[:,:,a,b]
            for j in rocc
                for i in rocc
                    temp = 0.0
                    temp += cache_iajb[i, j]
                    @views temp += dot(tiJaB_i[i,j,:,:], Wabef[a,b,:,:])
                    @views temp += dot(tiJaB_i[i,j,a,:], Fae[b,:])
                    @views temp += dot(tiJaB_i[j,i,b,:], Fae[a,:])
                    #for e in rvir
                    #    temp += tiJaB_i[i, j, a, e] * Fae[b, e]
                    #    temp += tiJaB_i[j, i, b, e] * Fae[a, e]
                    #    for f in rvir
                    #        temp += tiJaB_i[i, j, e, f] * Wabef[a, b, e, f]
                    #    end
                    #end
                    for m in rocc
                        temp -= tiJaB_i[i, m, a, b] * Fmi[m, j]
                        temp -= tiJaB_i[m, j, a, b] * Fmi[m, i]
                        for n in rocc
                            temp += tiJaB_i[m, n, a, b] * Wmnij[m, n, i, j]
                        end
                        for e in rvir
                            c1 = WmBeJ[m,b,e,j]
                            c2 = tiJaB_i[i,m,a,e]
                            c3 = WmBeJ[m,a,e,i]
                            c4 = tiJaB_i[j,m,b,e]
                            temp += (2 *( c2 * c1 + ( c4 * c3)) - (tiJaB_i[m, i, a, e] * c1) + (c2 * WmBEj[m, b, e, j])
                                     + (tiJaB_i[m, i, b, e] * WmBEj[m, a, e, j])
                                     + (tiJaB_i[m, j, a, e] * WmBEj[m, b, e, i])
                                     - (tiJaB_i[m, j, b, e] * c3)
                                     + (c4 * WmBEj[m, a, e, i]))
                        end
                    end
                    cache_tijab[i,j] = temp
                end
            end
            tiJaB_d[:,:,a,b] = cache_tijab ./ Dijab[:,:,a,b]
        end
    end
    end
    println("T2 in $dt")
    #tiJaB_i = nothing
    return convert(Array,tiJaB_d)
end
function form_Wmnij(minj, menf, tiJaB)
    dtt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Wmnij = zeros(dtt, nocc, nocc, nocc, nocc)
    form_Wmnij!(Wmnij, minj, menf, tiJaB)
    return Wmnij
end
function form_Wmnij!(Wmnij, minj, menf, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    for m in rocc
        for n in rocc
            cache_minj = minj[m,:,n,:]
            cache_menf = menf[m,:,n,:]
            for i in rocc
                for j in rocc
                    Wmnij[m, n, i, j] = cache_minj[i, j]
                    @views Wmnij[m, n, i, j] += dot(tiJaB[i,j,:,:],cache_menf[:,:])/2.0
                    #for f in rvir
                    #    @simd for e in rvir
                    #        Wmnij[m, n, i, j] += tiJaB[i, j, e, f] * cache_menf[e, f]/2.0 
                    #    end
                    #end
                end
            end
        end
    end
end

function form_Wabef(aebf, mnef, tiJaB)
    dt = eltype(aebf)
    nvir = size(tiJaB, 4)
    Wabef = zeros(dt, nvir, nvir, nvir, nvir)
    form_Wabef!(Wabef, aebf, mnef, tiJaB)
    return Wabef
end
function form_Wabef!(Wabef, aebf, menf, tiJaB)
    dtt = eltype(Wabef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    for f in rvir
        for e in rvir
            cache_aebf = aebf[:,e,:,f]
            cache_menf = menf[:,e,:,f]
            for b in rvir
                for a in rvir
                    #Wabef[a, b, e, f] = 0 
                    @views Wabef[a, b, e, f] = dot(tiJaB[:,:,a,b], cache_menf[:,:])/2.0
                    #for n in rocc
                    #    @simd for m in rocc
                    #        Wabef[a, b, e, f] += tiJaB[m, n, a, b] * cache_menf[m, n] 
                    #    end
                    #end
                    #Wabef[a, b, e, f] /= 2.0
                    Wabef[a, b, e, f] += cache_aebf[a, b]
                end
            end
        end
    end
end

function form_WmBeJ(mebj, iajb, tiJaB)
    dtt = eltype(iajb)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBeJ = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBeJ!(WmBeJ, mebj, iajb, tiJaB)
    return WmBeJ
end
function form_WmBeJ!(WmBeJ, mebj, iajb, tiJaB)
    dtt = eltype(WmBeJ)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    #WmBeJ .= 0 
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    for e in rvir
        for m in rocc
            cache_mebj = mebj[m,e,:,:]
            cache_iajb1 = iajb[m,e,:,:]
            cache_iajb2 = iajb[:,e,m,:]
            for b in rvir
                for j in rocc
                    WmBeJ[m, b, e, j] = 0#cache_mebj[b, j]
                    for f in rvir
                        @simd for n in rocc
                            WmBeJ[m, b, e, j] +=
                                cache_iajb1[n, f] *
                                (2 * tiJaB[n, j, f, b] - tiJaB[j, n, f, b])
                            WmBeJ[m, b, e, j] -= cache_iajb2[n, f] * tiJaB[n, j, f, b]
                        end
                    end
                    WmBeJ[m, b, e, j] /= 2.0
                    WmBeJ[m, b, e, j] += cache_mebj[b, j]
                end
            end
        end
    end
end


function form_WmBEj(nemf, mjbe, tiJaB)
    dtt = eltype(nemf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBEj = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBEj!(WmBEj, nemf, mjbe, tiJaB)
    return WmBEj
end
function form_WmBEj!(WmBEj, nemf, mjbe, tiJaB)
    dtt = eltype(WmBEj)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    for e in rvir
        for m in rocc
            cache_mjbe = mjbe[m,:,:,e]
            cache_nemf = nemf[:,e,m,:]
            for b in rvir
                for j in rocc
                    #WmBEj[m, b, e, j] = 0#-cache_mjbe[j, b]
                    @views WmBEj[m, b, e, j] = dot(tiJaB[j, :, :, b], cache_nemf[:, :])/2.0
                    #for f in rvir
                    #    @simd for n in rocc
                    #        WmBEj[m, b, e, j] += tiJaB[j, n, f, b] * cache_nemf[n, f] 
                    #    end
                    #end
                    #WmBEj[m, b, e, j] /= 2.0
                    WmBEj[m, b, e, j] += -cache_mjbe[j, b]
                end
            end
        end
    end
    return WmBEj
end
