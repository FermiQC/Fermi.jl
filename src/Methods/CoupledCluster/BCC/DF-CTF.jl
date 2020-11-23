using LinearAlgebra
using Fermi.DIIS

# Functions specific to DF-CTF
#
############ PREPROCESSING ################
"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function BCCD{T}(guess::BCCD{Tb},Alg::DFCTF) where { T <: AbstractFloat,
                                                    Tb <: AbstractFloat }
    refwfn = Fermi.HartreeFock.RHF()

    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    ints = refwfn.ints
    Fermi.Integrals.aux_ri!(ints)
    BCCD{T}(refwfn, guess, ints, Alg) 
end

function BCCD{T}(refwfn::RHF, guess::BCCD{Tb}, ints::IntegralHelper{Tc}, Alg::DFCTF) where { T <: AbstractFloat,
                                                                                               Tb <: AbstractFloat,
                                                                                              Tc <: AbstractFloat }
    d = [i - a for i = diag(ints["FOO"]), a = diag(ints["FVV"])]
    D = [i + j - a - b for i = diag(ints["FOO"]), j = diag(ints["FOO"]), a = diag(ints["FVV"]), b = diag(ints["FVV"])]
    newT1 = ints["FOV"]./d
    Bov = ints["BOV"]
    @tensor oovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
    newT2 = oovv ./ D
    BCCD{T}(refwfn, ints, newT1, newT2, Alg)
end
############ END PREPROCESSING ################
#
#
#
############# KERNEL FUNCTIONS ###################

function print_bcc_alg(Alg::DFCTF)
    @output "\n    • Computing BCC with the DF-CTF algorithm .\n\n"
end

function compute_integrals(ints,Alg::DFCTF)
    Fermi.Integrals.aux_ri!(ints)
    @output "Aux basis: {}\n" ints.bname["aux"]
    ints["BOV"]
    ints["BOO"]
    ints["BVV"]
end

function delete_integrals(ints,alg::DFCTF)
    delete!(ints.cache,"BOV")
    delete!(ints.cache,"BOO")
    delete!(ints.cache,"BVV")
    delete!(ints.cache,"F")
    delete!(ints.cache,"FOO")
    delete!(ints.cache,"FVV")
    delete!(ints.cache,"FOV")
end


function compute_oovv(ints,alg::DFCTF)
    Bov = ints["BOV"]
    @tensor oovv[i,j,a,b] := Bov[Q,i,a]*Bov[Q,j,b]
    oovv
end

function update_T1(T2::Array{T,4}, newT1::Array{T,2}, foo, fov, fvv, ints::IntegralHelper, alg::DFCTF) where { T <: AbstractFloat }
    newT1 .= 0
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
        newT1[i,a] += 2.0*fov[k,c]*T2[i,k,a,c]
        newT1[i,a] -= fov[k,c]*T2[k,i,a,c]
        newT1[i,a] -= T2[k,i,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += T2[l,k,a,c]*Vooov[k,l,i,c]
    end
end

function update_T2(T2::Array{T,4},newT2::Array{T,4},foo,fov,fvv,ints::IntegralHelper,alg::DFCTF) where T <: AbstractFloat
    newT2 .= 0
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
    T2_inter = zeros(T,dfsz,nvir,nvir)
    T2_slice = zeros(T,nvir,nvir)
    Tslice = zeros(T,nvir,nvir)
    for i=1:nocc, j=1:nocc
        T2_inter .= 0
        T2_slice .= 0
        @views Tslice .= T2[i,j,:,:]
        Fermi.contract!(T2_inter,Tslice,Bvv,"Qad","cd","Qca")
        Fermi.contract!(T2_slice,T2_inter,Bvv,"ab","Qad","Qdb")
        newT2[i,j,:,:] += T2_slice
    end

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>50x) begin
        newT2[i,j,a,b] += Voovv[i,j,a,b]
        newT2[i,j,a,b] += T2[k,l,a,b]*Voooo[i,j,k,l]
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
        P_OoVv[i,j,a,b] := -1.0*foo[i,k]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += fvv[c,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,i,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += 2.0*T2[i,k,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[i,k,a,c]*Vovov[j,c,k,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,j,a,c]*Vovov[i,c,k,b]
        P_OoVv[i,j,a,b] += T2[j,k,c,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T2[i,k,d,c]*T2[l,j,a,b]*Voovv[k,l,c,d]

        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end
end
#
#
############# END KERNEL FUNCTIONS ###################
