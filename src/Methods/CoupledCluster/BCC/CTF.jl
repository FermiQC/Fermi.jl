using LinearAlgebra
using Fermi.DIIS
# Functions specific to CTF
#
############ PREPROCESSING ################
"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function BCCD{T}(guess::BCCD{Tb},Alg::CTF) where { T <: AbstractFloat,
                                                    Tb <: AbstractFloat }
    refwfn = Fermi.HartreeFock.RHF()

    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    @output "Transforming Integrals..."
    tint = @elapsed begin
        refwfn.ints["OOOO"]
        refwfn.ints["OOOV"]
        refwfn.ints["OOVV"]
        refwfn.ints["OVOV"]
        refwfn.ints["OVVV"]
        refwfn.ints["VVVV"]
    end
    @output " done in {} s" tint
    BCCD{T}(refwfn, guess, Alg) 
end


function BCCD{T}(refwfn::RHF, guess::BCCD{Tb}, Alg::CTF) where { T <: AbstractFloat,
                                                                  Tb <: AbstractFloat }
    ints = refwfn.ints
    d = [i - a for i = diag(ints["FOO"]), a = diag(ints["FVV"])]
    D = [i + j - a - b for i = diag(ints["FOO"]), j = diag(ints["FOO"]), a = diag(ints["FVV"]), b = diag(ints["FVV"])]
    newT1 = ints["FOV"]./d
    newT2 = ints["OOVV"]./D
    o_small,v_small = size(guess.T1)
    o = 1:o_small
    v = 1:v_small
    newT1[o,v] .= Fermi.data(guess.T1)
    newT2[o,o,v,v] .= Fermi.data(guess.T2)
    BCCD{T}(refwfn, ints, newT1, newT2, Alg)
end
############ END PREPROCESSING ################
#
#
#
############# KERNEL FUNCTIONS ###################
function print_bcc_alg(alg::CTF)
    @output "Computing BCC with CTF algorithm.\n"
end

function compute_integrals(ints,alg::CTF)
    ints["OOOO"]
    ints["OOOV"]
    ints["OOVV"]
    ints["OVOV"]
    ints["OVVV"]
    ints["VVVV"]
end

function delete_integrals(ints,alg::CTF)
    delete!(ints.cache,"OOOO")
    delete!(ints.cache,"OOOV")
    delete!(ints.cache,"OOVV")
    delete!(ints.cache,"OVOV")
    delete!(ints.cache,"OVVV")
    delete!(ints.cache,"VVVV")
    delete!(ints.cache,"FOO")
    delete!(ints.cache,"FVV")
    delete!(ints.cache,"FOV")
end


function compute_oovv(ints,alg::CTF)
    ints["OOVV"]
end

function update_T1(T2::Array{T,4}, newT1::Array{T,2}, foo, fov, fvv, ints::IntegralHelper, alg::CTF) where T <: AbstractFloat
    fill!(newT1, 0.0)
    Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = ints["OOOO"], ints["OOOV"], ints["OOVV"], ints["OVOV"], ints["OVVVV"], ints["VVVV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
        newT1[i,a] += fov[i,a]
        newT1[i,a] += 2.0*fov[k,c]*T2[i,k,a,c]
        newT1[i,a] -= fov[k,c]*T2[k,i,a,c]
        newT1[i,a] -= T2[k,i,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += T2[l,k,a,c]*Vooov[k,l,i,c]
    end
end

function update_T2(T2::Array{T,4},newT2::Array{T,4},foo,fov,fvv,ints::IntegralHelper, alg::CTF) where T <: AbstractFloat
    fill!(newT2, 0.0)
    Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = ints["OOOO"], ints["OOOV"], ints["OOVV"], ints["OVOV"], ints["OVVVV"], ints["VVVV"]

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
        newT2[i,j,a,b] += Voovv[i,j,a,b]
        newT2[i,j,a,b] += T2[i,j,c,d]*Vvvvv[c,d,a,b]
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
