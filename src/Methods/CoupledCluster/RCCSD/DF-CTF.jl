using LinearAlgebra
using Fermi.DIIS

# Functions specific to DF-CTF
#
############ PREPROCESSING ################
"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function RCCSD{T}(guess::RCCSD{Tb},Alg::DFCTF) where { T <: AbstractFloat,
                                                    Tb <: AbstractFloat }
    refwfn = Fermi.HartreeFock.RHF()

    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    ints = refwfn.ints
    Fermi.Integrals.aux_ri!(ints)
    RCCSD{T}(refwfn, guess, ints, Alg) 
end

function RCCSD{T}(refwfn::RHF, guess::RCCSD{Tb}, ints::IntegralHelper{Tc}, Alg::DFCTF) where { T <: AbstractFloat,
                                                                                               Tb <: AbstractFloat,
                                                                                              Tc <: AbstractFloat }
    d = [i - a for i = diag(ints["FOO"]), a = diag(ints["FVV"])]
    D = [i + j - a - b for i = diag(ints["FOO"]), j = diag(ints["FOO"]), a = diag(ints["FVV"]), b = diag(ints["FVV"])]
    newT1 = ints["FOV"]./d
    Bov = ints["BOV"]
    @tensor oovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
    newT2 = oovv ./ D
    RCCSD{T}(refwfn, ints, newT1, newT2, Alg)
end
############ END PREPROCESSING ################
#
#
#
############# KERNEL FUNCTIONS ###################

function print_alg(Alg::DFCTF)
    @output "\n    • Computing CCSD with the DF-CTF algorithm .\n\n"
end

function compute_integrals(ints,Alg::DFCTF)
    Fermi.Integrals.aux_ri!(ints)
    @output "Aux basis: {}\n" ints.bname["aux"]
    ints["BOV"]
    ints["BOO"]
    ints["BVV"]
end

function compute_oovv(ints,alg::DFCTF)
    Bov = ints["BOV"]
    @tensor oovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
    oovv
end

function update_T1(T1::Array{T,2}, T2::Array{T,4}, newT1::Array{T,2}, foo, fov, fvv, ints::IntegralHelper, alg::DFCTF) where { T <: AbstractFloat }
    Bov = ints["BOV"]
    Boo = ints["BOO"]
    Bvv = ints["BVV"]
    @tensor begin
        Voooo[i,j,k,l] := Boo[Q,i,k]*Boo[Q,j,l]
        Vooov[i,j,k,a] := Boo[Q,i,k]*Bov[Q,j,a]
        Voovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
        Vovov[i,a,j,b] := Boo[Q,i,j]*Bvv[Q,a,b]
    end
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x,Q=>50x) begin
        newT1[i,a] += fov[i,a]
        newT1[i,a] -= foo[i,k]*T1[k,a]
        newT1[i,a] += fvv[c,a]*T1[i,c]
        newT1[i,a] -= fov[k,c]*T1[i,c]*T1[k,a]
        newT1[i,a] += 2.0*fov[k,c]*T2[i,k,a,c]
        newT1[i,a] -= fov[k,c]*T2[k,i,a,c]
        newT1[i,a] -= T1[k,c]*Vovov[i,c,k,a]
        newT1[i,a] += 2.0*T1[k,c]*Voovv[k,i,c,a]
        newT1[i,a] -= T2[k,i,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += T2[l,k,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += -2.0*T1[k,c]*T1[l,a]*Vooov[l,k,i,c]
        newT1[i,a] -= T1[k,c]*T1[i,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Bov[Q,k,c]*Bvv[Q,a,d]
        newT1[i,a] += T1[k,c]*T1[l,a]*Vooov[k,l,i,c]
        newT1[i,a] += -2.0*T1[k,c]*T2[i,l,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[k,c]*T2[l,i,a,d]*Voovv[k,l,c,d]
        newT1[i,a] += T1[k,c]*T2[l,i,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[i,c]*T2[l,k,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += T1[i,c]*T2[l,k,a,d]*Voovv[k,l,c,d]
        newT1[i,a] += -2.0*T1[l,a]*T2[i,k,d,c]*Voovv[k,l,c,d]
        newT1[i,a] += T1[l,a]*T2[i,k,c,d]*Voovv[k,l,c,d]
        newT1[i,a] += T1[k,c]*T1[i,d]*T1[l,a]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[k,c]*T1[i,d]*T1[l,a]*Voovv[k,l,c,d]
        newT1[i,a] += 4.0*T1[k,c]*T2[i,l,a,d]*Voovv[k,l,c,d]
    end
end

function update_T2(T1::Array{T,2},T2::Array{T,4},newT2::Array{T,4},foo,fov,fvv,ints::IntegralHelper,alg::DFCTF) where T <: AbstractFloat
    Bov = ints["BOV"]
    Boo = ints["BOO"]
    Bvv = ints["BVV"]
    @tensor begin
        Voooo[i,j,k,l] := Boo[Q,i,k]*Boo[Q,j,l]
        Vooov[i,j,k,a] := Boo[Q,i,k]*Bov[Q,j,a]
        Voovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
        Vovov[i,a,j,b] := Boo[Q,i,j]*Bvv[Q,a,b]
    end
    dfsz,nocc,nvir = size(Bov)
    @tensor _T[i,j,c,d] := T1[i,c]*T1[j,d] + T2[i,j,c,d]
    T2_inter = zeros(T,dfsz,nvir,nvir)
    T2_slice = zeros(T,nvir,nvir)
    Tslice = zeros(T,nvir,nvir)
    for i=1:nocc, j=1:nocc
        T2_inter .= 0
        T2_slice .= 0
        @views Tslice .= _T[i,j,:,:]
        Fermi.contract!(T2_inter,Tslice,Bvv,"Qad","cd","Qca")
        Fermi.contract!(T2_slice,T2_inter,Bvv,"ab","Qad","Qdb")
        newT2[i,j,:,:] += T2_slice
    end

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>50x) begin
        newT2[i,j,a,b] += Voovv[i,j,a,b]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] += T2[k,l,a,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,a]*Bov[Q,k,c]*Bvv[Q,b,d]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,b]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT2[i,j,a,b] += T1[i,c]*T1[k,a]*T1[l,b]*Vooov[l,k,j,c]
        newT2[i,j,a,b] += T1[j,c]*T1[k,a]*T1[l,b]*Vooov[k,l,i,c]
        newT2[i,j,a,b] += T2[k,l,a,c]*T2[i,j,d,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[l,j,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[l,k,a,c]*T2[i,j,d,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,d,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T2[i,k,a,c]*T2[l,j,b,d]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[j,l,b,d]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[k,i,a,c]*T2[j,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[i,j,a,c]*T2[l,k,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,j,a,c]*T2[k,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[k,j,a,c]*T2[i,l,d,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += 4.0*T2[i,k,a,c]*T2[j,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[i,j,d,c]*T2[l,k,a,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T1[k,a]*T1[l,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T2[l,k,a,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*T2[i,j,d,c]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] := -1.0*foo[i,k]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += fvv[c,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[k,b]*Vooov[j,i,k,a]
        P_OoVv[i,j,a,b] += T1[j,c]*Bov[Q,i,a]*Bvv[Q,c,b]
        P_OoVv[i,j,a,b] += -1.0*fov[k,c]*T1[i,c]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += -1.0*fov[k,c]*T1[k,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,i,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,a]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,b]*Vovov[j,c,k,a]
        P_OoVv[i,j,a,b] += 2.0*T2[i,k,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[i,k,a,c]*Vovov[j,c,k,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,j,a,c]*Vovov[i,c,k,b]
        P_OoVv[i,j,a,b] += -2.0*T1[l,b]*T2[i,k,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[k,i,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,d,b]*Bov[Q,k,c]*Bvv[Q,a,d]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[k,i,a,d]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,a,d]*Bov[Q,k,c]*Bvv[Q,b,d]
        P_OoVv[i,j,a,b] += 2.0*T1[k,c]*T2[i,j,a,d]*Bov[Q,k,c]*Bvv[Q,b,d]
        P_OoVv[i,j,a,b] += T1[j,c]*T2[l,k,a,b]*Vooov[l,k,i,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[i,k,a,c]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,a]*T2[i,j,d,c]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += T1[k,a]*T2[i,l,c,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += 2.0*T1[j,c]*T2[i,k,a,d]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] -= 1.0*T1[k,c]*T2[i,j,a,d]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += T1[k,c]*T2[i,l,a,b]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T2[i,l,a,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += T2[j,k,c,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[i,c]*T1[k,a]*T2[j,l,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[l,b]*T2[k,j,a,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T2[i,k,d,c]*T2[l,j,a,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end
end
#
#
############# END KERNEL FUNCTIONS ###################
