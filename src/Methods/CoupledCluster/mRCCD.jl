"""
    Fermi.CoupledCluster.mRCCD

performs coupled cluster doubles (CCD) computations using explicit matrix multiplications.

## Functions

    Fermi.CoupledCluster.RCCD.do_rccd
"""
module mRCCD
using Base.Threads
using Fermi.Wavefunction
#using CuArrays
using Fermi
using Fermi.Output
using TensorOperations
using LinearAlgebra
using Fermi.Transformation
using Fermi.IntegralTransformation
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
    do_diis = true
    #DIIS ripped straight from psi4numpy
    Fermi.CoupledCluster.print_header()
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    epsa = refWfn.epsa
    T = eltype(refWfn.ao_eri)
    oovv,ovov,ovvo,oooo,vvvv = make_rccd_integrals(refWfn) 
    T2 = zeros(T, nocc, nocc, nvir, nvir)
    Dijab = form_Dijab(T2, epsa)
    T2_init!(T2, oovv, Dijab)
    if doprint
        println("@RMP2 ", ccenergy(T2, oovv))
    end
    oovv = Array{T}(oovv)
    ovov = Array{T}(ovov)
    ovvo = Array{T}(ovvo)
    oooo = Array{T}(oooo)
    vvvv = Array{T}(vvvv)
    Fae = form_Fae(T2, oovv)
    Fmi = form_Fmi(T2, oovv)
    Wmnij = form_Wmnij(oooo, oovv, T2)
    Wabef = form_Wabef(vvvv, oovv, T2)
    WmBeJ = form_WmBeJ(ovvo, oovv, T2)
    WmBEj = form_WmBEj(oovv, ovov, T2)
    diis_vals_t2 = [convert(Array{Float32},deepcopy(T2))]
    diis_errors = [zeros(Float32,0)]
    diis_size = 0
    max_diis = 6
    dt = @elapsed for i in 0:maxit-1 #TODO: implement RMS check
        if do_diis
            T2,rms = cciter(
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
                diis_vals_t2,
                diis_errors,
                max_diis,
                i
               )
        else
            T2,rms = cciter(
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
                WmBEj
            )
        end
        if rms < 1E-7 #TODO have proper kwargs
            break
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

function make_rccd_integrals(wfn)
    oovv = get_eri(wfn,"OOVV")
    vvvv = get_eri(wfn,"VVVV")
    ovvo = get_eri(wfn,"OVVO")
    ovov = get_eri(wfn,"OVOV")
    oooo = get_eri(wfn,"OOOO")
    return oovv,ovov,ovvo,oooo,vvvv
end

function ccenergy(tiJaB, ijab)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    @tensor ecc = ijab[i,j,a,b]*2*tiJaB[i,j,a,b] - ijab[i,j,a,b]*tiJaB[j,i,a,b]
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
    diis_vals_t2,
    diis_errors,
    max_diis,
    i
)
    o = size(tiJaB_i,1)
    v = size(tiJaB_i,3)
    form_Wabef!(Wabef, vvvv, oovv, tiJaB_i)
    form_WmBeJ!(WmBeJ, ovvo, oovv, tiJaB_i)
    form_WmBEj!(WmBEj, oovv, ovov, tiJaB_i)
    form_Wmnij!(Wmnij, oooo, oovv, tiJaB_i)
    form_Fae!(Fae, tiJaB_i, oovv)
    form_Fmi!(Fmi, tiJaB_i, oovv)

    tiJaB_d = form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, oovv, Dijab)
    e = ccenergy(tiJaB_d,oovv)
    #println("$e")
    @output "@CCD {:20.17f}" ccenergy(tiJaB_d,oovv)
    e = ccenergy(tiJaB_d,oovv)
    push!(diis_vals_t2,convert(Array{Float32},deepcopy(tiJaB_d)))
    error_t2 = reshape(tiJaB_d - tiJaB_i,o^2*v^2)
    #error_t2 = flat(tiJaB_d - tiJaB_i)
    #println(norm(error_t2))
    push!(diis_errors,convert(Array{Float32},deepcopy(error_t2)))
    i == 0 ? deleteat!(diis_errors,1) : nothing

    diis_size = length(diis_errors)
    if i >= 0
        if length(diis_vals_t2) > max_diis
            deleteat!(diis_vals_t2,1)
            deleteat!(diis_errors,1)
        end
        diis_size = length(diis_vals_t2) - 1
        B = ones(Float32,diis_size+1, diis_size+1)* -1
        B[end,end] = 0
        
        for (n1, e1) in enumerate(diis_errors[1:end])
            for (n2, e2) in enumerate(diis_errors[1:end])
                B[n1,n2] = dot(e1,e2)
            end
        end
        E = size(B,1)
        B[1:E-1,1:E-1] ./= maximum(abs.(B[1:E-1,1:E-1]))
        resid = zeros(Float32,diis_size+1)
        resid[end] = -1
        #LAPACK.gesv!(B,resid)
        ci = inv(B)*resid
        tiJaB_d .= 0
        for num in 1:diis_size
            tiJaB_d .+= convert(Float32,ci[num])*diis_vals_t2[num+1]
        end
    end


    r2 = norm(error_t2)#sqrt(sum((tiJaB_d - tiJaB_i).^2))/length(tiJaB_d)
    return tiJaB_d, r2
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
    e = ccenergy(tiJaB_d,oovv)
    #println("@CCD $e")
    r2 = sqrt(sum((tiJaB_d - tiJaB_i).^2))/length(tiJaB_d)
    return tiJaB_d,r2
end

function T2_init!(tiJaB, iajb, Dijab)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    #tiJaB .= permutedims(iajb,(1,3,2,4)) ./ Dijab
    tiJaB .= zeros(eltype(tiJaB),size(iajb)) ./ Dijab
    #for i in rocc
    #    for j in rocc
    #        cache = iajb[i,:,j,:]
    #        for a in rvir
    #            for b in rvir
    #                tiJaB[i,j,a,b] = cache[a,b] / Dijab[i,j,a,b]
    #            end
    #        end
    #    end
    #end
end

function form_Fae(tiJaB, menf)
    dt = eltype(menf)
    nvir = size(tiJaB, 4)
    Fae = zeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, menf)
    return Fae
end
function form_Fae!(Fae, tiJaB, mnef)
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _tiJaB = 2*reshape(permutedims(tiJaB,(3,1,2,4)),v,o^2*v)
    _tiJaB -= reshape(permutedims(tiJaB,(3,2,1,4)),v,o^2*v)
    _mnef = reshape(permutedims(mnef,(1,2,4,3)),o^2*v,v)
    T = eltype(tiJaB)
    c1::T = -1.0
    c2::T = 0.0
    BLAS.gemm!('N','N',c1,_tiJaB,_mnef,c2,Fae)
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
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    #@cast _tiJaB[(n,e,f),(i)] := 2*tiJaB[i,n,e,f] - tiJaB[i,n,f,e]
    _tiJaB = 2*reshape(permutedims(tiJaB,(2,3,4,1)),o*v^2,o)
    _tiJaB -= reshape(permutedims(tiJaB,(2,4,3,1)),o*v^2,o)
    #@cast _mnef[(m),(n,e,f)] := mnef[m,n,e,f]
    _mnef = reshape(mnef,o,o*v^2)
    #mul!(Fmi,_mnef,_tiJaB)
    T = eltype(tiJaB)
    c1::T = 1.0
    c2::T = 0.0
    BLAS.gemm!('N','N',c1,_mnef,_tiJaB,c2,Fmi)
    return Fmi
end

@fastmath @inbounds function form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, ijab, Dijab)
    chonk = 1
    T = eltype(tiJaB_i)
    _one::T = 1.0
    _p5::T = 0.5
    _zero::T = 0.0
    _neg::T = -1.0
    tiJaB_d = zeros(T,size(ijab[:,:,:,:]))
    tiJaB_d .= ijab
    o = size(ijab,1)
    v = size(ijab,3)
    _tiJaB_i = reshape(tiJaB_i,o^2*v,v)
    _tiJaB_d = reshape(tiJaB_d,o^2*v,v)
    BLAS.gemm!('N','T',_one,_tiJaB_i,Fae,_one,_tiJaB_d)
    tiJaB_d .= reshape(_tiJaB_d,o,o,v,v)

    @tensor _tiJaB_d[j,i,b,a] := tiJaB_d[i,j,a,b]
    #_tiJaB_d = reshape(permutedims(tiJaB_d,(2,1,4,3)),o^2*v,v)
    _tiJaB_d = reshape(_tiJaB_d,o^2*v,v)
    _tiJaB_i = reshape(tiJaB_i,o^2*v,v)
    BLAS.gemm!('N','T',_one,_tiJaB_i,Fae,_one,_tiJaB_d)
    tiJaB_d .= permutedims(reshape(_tiJaB_d,o,o,v,v),(2,1,4,3))

    @tensor _tiJaB_d[1,3,4,2] := tiJaB_d[1,2,3,4]
    _tiJaB_d = reshape(_tiJaB_d,o*v^2,o)
    _tiJaB_i = reshape(permutedims(tiJaB_i,(2,1,3,4)),o,o*v^2)
    BLAS.gemm!('T','N',_neg*_one,_tiJaB_i,Fmi,_one,_tiJaB_d)
    tiJaB_d .= permutedims(reshape(transpose(_tiJaB_d),o,o,v,v),(2,1,3,4))

    #_tiJaB_d = reshape(permutedims(tiJaB_d,(2,3,4,1)),o*v^2,o)
    @tensor _tiJaB_d2[2,3,4,1] := tiJaB_d[1,2,3,4]
    _tiJaB_d2 = reshape(_tiJaB_d2,o*v^2,o)
    _tiJaB_i2 = reshape(tiJaB_i,o,o*v^2)
    BLAS.gemm!('T','N',_neg*_one,_tiJaB_i2,Fmi,_one,_tiJaB_d2)
    tiJaB_d .= reshape(transpose(_tiJaB_d2),o,o,v,v)

    _tiJaB_d = reshape(tiJaB_d,o^2,v^2)
    _tiJaB_i = reshape(tiJaB_i,o^2,v^2)
    _Wmnij = reshape(Wmnij,o^2,o^2)
    BLAS.gemm!('T','N',_one,_Wmnij,_tiJaB_i,_one,_tiJaB_d)

    _Wabef = reshape(Wabef,v^2,v^2)
    BLAS.gemm!('N','T',_one,_tiJaB_i,_Wabef,_one,_tiJaB_d)
    #GSGEMM!(_tiJaB_d,1.0f0,_tiJaB_i,_Wabef,1.0f0,size(_tiJaB_i,2),size(_tiJaB_i,2))
    tiJaB_d = reshape(_tiJaB_d,o,o,v,v)

    scr1 = Array{T}(undef,o*v,o*v)
    scr2 = Array{T}(undef,o*v,o*v)
    scr3 = Array{T}(undef,o*v,o*v)

    @tensor _tiJaB_d[1,3,4,2] := tiJaB_d[1,2,3,4]
    _tiJaB_d = reshape(_tiJaB_d,o*v,v*o)
    @tensor _scr2[1,3,2,4] := tiJaB_i[1,2,3,4]
    scr2 = reshape(_scr2,o*v,o*v)
    @tensor _scr3[1,3,2,4] := WmBeJ[1,2,3,4]
    scr3 = 2*reshape(_scr3,o*v,v*o)

    @tensor int[1,3,2,4] := WmBEj[1,2,3,4]
    int = reshape(int,o*v,v*o)
    scr3 += int
    BLAS.gemm!('N','N',_one,scr2,scr3,_one,_tiJaB_d)

    scr3 -= int
    scr3 *= 0.5f0
    scr2 = reshape(permutedims(tiJaB_i,(1,4,2,3)),o*v,o*v)
    BLAS.gemm!('T','N',_neg*_one,scr2,scr3,_one,_tiJaB_d)

    tiJaB_d = permutedims(reshape(_tiJaB_d,o,v,v,o),(1,4,2,3))

    _tiJaB_d = reshape(permutedims(tiJaB_d,(1,4,3,2)),o*v,v*o)
    scr2 = reshape(permutedims(tiJaB_i,(1,4,2,3)),o*v,o*v)
    scr4 = reshape(permutedims(WmBEj,(1,3,2,4)),o*v,o*v)
    tmp = Array{T}(undef,size(_tiJaB_d))
    BLAS.gemm!('T','N',_one,scr2,scr4,_zero,tmp)
    _tiJaB_d += tmp
    tiJaB_d = permutedims(reshape(_tiJaB_d,o,v,v,o),(1,4,3,2))
    

    scr1 = reshape(permutedims(tiJaB_d,(2,3,4,1)),o*v,v*o)
    scr1 += tmp
    tiJaB_d = permutedims(reshape(scr1,o,v,v,o),(4,1,2,3))

    scr1 = reshape(permutedims(tiJaB_d,(2,4,3,1)),o*v,v*o)
    scr2 = reshape(permutedims(tiJaB_i,(1,3,2,4)),o*v,o*v)
    BLAS.gemm!('N','N',_one,scr2,scr4,_one,scr1)

    scr2 *= 2
    scr2 -= reshape(permutedims(tiJaB_i,(2,3,1,4)),o*v,o*v)
    BLAS.gemm!('N','N',_one,scr2,scr3,_one,scr1)

    tiJaB_d = permutedims(reshape(scr1,o,v,v,o),(4,1,3,2)) ./ Dijab

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
    Wmnij .= mnij
    T = eltype(tiJaB)
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _Wmnij = reshape(Wmnij,o^2,o^2)
    _tiJaB = reshape(tiJaB,o^2,v^2)
    _mnef = reshape(mnef,o^2,v^2)
    c1::T = 0.5
    c2::T = 1.0
    BLAS.gemm!('N','T',c1,_mnef,_tiJaB,c2,_Wmnij)
    Wmnij .= reshape(_Wmnij,o,o,o,o)
end

function form_Wabef(aebf, mnef, tiJaB)
    dt = eltype(aebf)
    nvir = size(tiJaB, 4)
    Wabef = zeros(dt, nvir, nvir, nvir, nvir)
    form_Wabef!(Wabef, aebf, mnef, tiJaB)
    return Wabef
end
function form_Wabef!(Wabef, abef, mnef, tiJaB)
    Wabef .= abef
    T = eltype(tiJaB)
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _Wabef = reshape(Wabef,v^2,v^2)
    _tiJaB = reshape(tiJaB,o^2,v^2)
    _mnef = reshape(mnef,o^2,v^2)
    c1::T = 0.5
    c2::T = 1.0
    BLAS.gemm!('T','N',c1,_tiJaB,_mnef,c2,_Wabef)
    Wabef = reshape(_Wabef,v,v,v,v)
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
    WmBeJ .= mbej
    T = eltype(tiJaB)
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _WmBeJ = reshape(permutedims(WmBeJ,(1,3,4,2)),o*v,o*v)
    _tiJaB1 = 2*reshape(permutedims(tiJaB,(1,3,2,4)),o*v,o*v)
    _tiJaB1 -= reshape(permutedims(tiJaB,(2,3,1,4)),o*v,o*v)
    _mnef1 = reshape(permutedims(mnef,(1,3,2,4)),o*v,o*v)

    _mnef2 = reshape(permutedims(mnef,(1,4,2,3)),o*v,o*v)
    _tiJaB2 = reshape(permutedims(tiJaB,(1,3,2,4)),o*v,o*v)
    _neg::T = -1.0
    c1::T = 0.5
    c2::T = 1.0
    BLAS.gemm!('N','N',c1,_mnef1,_tiJaB1,c2,_WmBeJ)
    BLAS.gemm!('T','N',_neg*c1,_mnef2,_tiJaB2,c2,_WmBeJ)
    WmBeJ .= permutedims(reshape(_WmBeJ,o,v,o,v),(1,4,2,3))
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
    WmBEj .= -permutedims(mbje,(1,2,4,3))
    T = eltype(tiJaB)
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _tiJaB = reshape(permutedims(tiJaB,(2,3,1,4)),o*v,o*v)
    _nmef = reshape(permutedims(nmef,(2,3,1,4)),o*v,o*v)
    _WmBEj = (reshape(permutedims(WmBEj,(1,3,4,2)),o*v,o*v))
    c1::T = 0.5
    c2::T = 1.0
    BLAS.gemm!('N','N',c1,_nmef,_tiJaB,c2,_WmBEj)
    WmBEj .= permutedims(reshape(_WmBEj,o,v,o,v),(1,4,2,3))
    return WmBEj
end
end #module
