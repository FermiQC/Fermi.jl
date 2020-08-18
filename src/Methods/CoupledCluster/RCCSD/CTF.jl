using LinearAlgebra
using Fermi.DIIS
# Functions specific to CTF
#
############ PREPROCESSING ################
"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function RCCSD{T}(guess::RCCSD{Tb},Alg::CTF) where { T <: AbstractFloat,
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
    RCCSD{T}(refwfn, guess, Alg) 
end


function RCCSD{T}(refwfn::RHF, guess::RCCSD{Tb}, Alg::CTF) where { T <: AbstractFloat,
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
    RCCSD{T}(refwfn, ints, newT1, newT2, Alg)
end
############ END PREPROCESSING ################
#
#
#
############# KERNEL FUNCTIONS ###################
function print_alg(alg::CTF)
    @output "Computing RCCSD with CTF algorithm.\n"
end

function compute_integrals(ints,alg::CTF)
    ints["OOOO"]
    ints["OOOV"]
    ints["OOVV"]
    ints["OVOV"]
    ints["OVVV"]
    ints["VVVV"]
end
function update_energy(T1::Array{T, 2}, T2::Array{T, 4}, f::Array{T,2}, ints, alg::CTF) where { T <: AbstractFloat }

    Voovv = ints["OOVV"]
    @tensoropt (k=>x, l=>x, c=>100x, d=>100x)  begin
        CC_energy = 2.0*f[k,c]*T1[k,c]
        B[l,c,k,d] := -1.0*T1[l,c]*T1[k,d]
        B[l,c,k,d] += -1.0*T2[l,k,c,d]
        B[l,c,k,d] += 2.0*T2[k,l,c,d]
        CC_energy += B[l,c,k,d]*Voovv[k,l,c,d]
        CC_energy += 2.0*T1[l,c]*T1[k,d]*Voovv[l,k,c,d]
    end
    
    return CC_energy
end

function compute_oovv(ints,alg::CTF)
    ints["OOVV"]
end

function update_T1(T1::Array{T,2}, T2::Array{T,4}, newT1::Array{T,2}, foo, fov, fvv, ints::IntegralHelper, alg::CTF) where T <: AbstractFloat
    fill!(newT1, 0.0)
    Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = ints["OOOO"], ints["OOOV"], ints["OOVV"], ints["OVOV"], ints["OVVVV"], ints["VVVV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
        newT1[i,a] += fov[i,a]
        newT1[i,a] -= foo[i,k]*T1[k,a]
        newT1[i,a] += fvv[c,a]*T1[i,c]
        newT1[i,a] -= fov[k,c]*T1[i,c]*T1[k,a]
        newT1[i,a] += 2.0*fov[k,c]*T2[i,k,a,c]
        newT1[i,a] -= fov[k,c]*T2[k,i,a,c]
        newT1[i,a] -= T1[k,c]*Vovov[i,c,k,a]
        newT1[i,a] += 2.0*T1[k,c]*Voovv[k,i,c,a]
        newT1[i,a] -= T2[k,i,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += T2[l,k,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += -2.0*T1[k,c]*T1[l,a]*Vooov[l,k,i,c]
        newT1[i,a] -= T1[k,c]*T1[i,d]*Vovvv[k,a,d,c]
        newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Vovvv[k,a,c,d]
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

function update_T2(T1::Array{T,2},T2::Array{T,4},newT2::Array{T,4},foo,fov,fvv,ints::IntegralHelper, alg::CTF) where T <: AbstractFloat
    fill!(newT2, 0.0)
    Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = ints["OOOO"], ints["OOOV"], ints["OOVV"], ints["OVOV"], ints["OVVVV"], ints["VVVV"]

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
        newT2[i,j,a,b] += Voovv[i,j,a,b]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*Vvvvv[c,d,a,b]
        newT2[i,j,a,b] += T2[i,j,c,d]*Vvvvv[c,d,a,b]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] += T2[k,l,a,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,a]*Vovvv[k,b,c,d]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,b]*Vovvv[k,a,d,c]
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
        P_OoVv[i,j,a,b] += T1[j,c]*Vovvv[i,c,a,b]
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
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,d,b]*Vovvv[k,a,c,d]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[k,i,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,a,d]*Vovvv[k,b,c,d]
        P_OoVv[i,j,a,b] += T1[j,c]*T2[l,k,a,b]*Vooov[l,k,i,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[i,k,a,c]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,a]*T2[i,j,d,c]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += T1[k,a]*T2[i,l,c,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += 2.0*T1[j,c]*T2[i,k,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,c]*T2[i,j,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += 2.0*T1[k,c]*T2[i,j,a,d]*Vovvv[k,b,c,d]
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
