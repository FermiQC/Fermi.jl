function cc_update_energy(T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,E,O}, 
                          alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Vovov = moints["OVOV"]
    @tensoropt (k=>x, l=>x, c=>100x, d=>100x)  begin
        B[l,c,k,d] := 2.0*T2[k,l,c,d]
        B[l,c,k,d] -= T1[l,c]*T1[k,d]
        B[l,c,k,d] -= T2[l,k,c,d]
        CC_energy = B[l,c,k,d]*Vovov[k,c,l,d]
        CC_energy += 2.0*T1[l,c]*T1[k,d]*Vovov[l,c,k,d]
    end
    
    return CC_energy
end

function cc_update_T1!(newT1::AbstractArray{T,2}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,E,O}, 
                      alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Vooov, Voovv, Vovov, Vovvv = moints["OOOV"], moints["OOVV"], moints["OVOV"], moints["OVVV"]

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
        newT1[i,a] -= T1[k,c]*Voovv[i,k,c,a]
        newT1[i,a] += 2.0*T1[k,c]*Vovov[k,c,i,a]
        newT1[i,a] -= T2[k,i,c,d]*Vovvv[k,d,a,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Vovvv[k,d,a,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,i,l,c]
        newT1[i,a] += T2[l,k,a,c]*Vooov[k,i,l,c]
        newT1[i,a] += -2.0*T1[k,c]*T1[l,a]*Vooov[l,i,k,c]
        newT1[i,a] -= T1[k,c]*T1[i,d]*Vovvv[k,d,a,c]
        newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Vovvv[k,c,a,d]
        newT1[i,a] += T1[k,c]*T1[l,a]*Vooov[k,i,l,c]
        newT1[i,a] += -2.0*T1[k,c]*T2[i,l,a,d]*Vovov[l,c,k,d]
        newT1[i,a] += -2.0*T1[k,c]*T2[l,i,a,d]*Vovov[k,c,l,d]
        newT1[i,a] += T1[k,c]*T2[l,i,a,d]*Vovov[l,c,k,d]
        newT1[i,a] += -2.0*T1[i,c]*T2[l,k,a,d]*Vovov[l,c,k,d]
        newT1[i,a] += T1[i,c]*T2[l,k,a,d]*Vovov[k,c,l,d]
        newT1[i,a] += -2.0*T1[l,a]*T2[i,k,d,c]*Vovov[k,c,l,d]
        newT1[i,a] += T1[l,a]*T2[i,k,c,d]*Vovov[k,c,l,d]
        newT1[i,a] += T1[k,c]*T1[i,d]*T1[l,a]*Vovov[l,c,k,d]
        newT1[i,a] += -2.0*T1[k,c]*T1[i,d]*T1[l,a]*Vovov[k,c,l,d]
        newT1[i,a] += 4.0*T1[k,c]*T2[i,l,a,d]*Vovov[k,c,l,d]
    end
end

function cc_update_T2!(newT2::AbstractArray{T,4}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,E,O}, 
                    alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Bvv = moints["BVV"]
    Voooo, Vooov, Voovv, Vovov, Vovvv = moints["OOOO"], moints["OOOV"], moints["OOVV"], moints["OVOV"], moints["OVVV"]

    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x, Q=>50x) begin
        newT2[i,j,a,b] += Vovov[i,a,j,b]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*Bvv[Q,c,a]*Bvv[Q,d,b]
        newT2[i,j,a,b] += T2[i,j,c,d]*Bvv[Q,c,a]*Bvv[Q,d,b]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*Voooo[i,k,j,l]
        newT2[i,j,a,b] += T2[k,l,a,b]*Voooo[i,k,j,l]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,a]*Vovvv[k,c,b,d]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,b]*Vovvv[k,d,a,c]
        newT2[i,j,a,b] += T1[i,c]*T1[k,a]*T1[l,b]*Vooov[l,j,k,c]
        newT2[i,j,a,b] += T1[j,c]*T1[k,a]*T1[l,b]*Vooov[k,i,l,c]
        newT2[i,j,a,b] += T2[k,l,a,c]*T2[i,j,d,b]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[l,j,b,d]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += -2.0*T2[l,k,a,c]*T2[i,j,d,b]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,d,b]*Vovov[l,c,k,d]
        newT2[i,j,a,b] += T2[i,k,a,c]*T2[l,j,b,d]*Vovov[l,c,k,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[j,l,b,d]*Vovov[l,c,k,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,b,d]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += -2.0*T2[k,i,a,c]*T2[j,l,b,d]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += T2[i,j,a,c]*T2[l,k,b,d]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += -2.0*T2[i,j,a,c]*T2[k,l,b,d]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += T2[k,j,a,c]*T2[i,l,d,b]*Vovov[l,c,k,d]
        newT2[i,j,a,b] += 4.0*T2[i,k,a,c]*T2[j,l,b,d]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += T2[i,j,d,c]*T2[l,k,a,b]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T1[k,a]*T1[l,b]*Vovov[k,c,l,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T2[l,k,a,b]*Vovov[l,c,k,d]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*T2[i,j,d,c]*Vovov[l,c,k,d]
        P_OoVv[i,j,a,b] := T1[j,c]*Vovvv[i,a,c,b]
        P_OoVv[i,j,a,b] -= T1[k,b]*Vooov[j,k,i,a]
        P_OoVv[i,j,a,b] -= T2[k,i,a,c]*Vovov[k,c,j,b]
        P_OoVv[i,j,a,b] -= T1[i,c]*T1[k,a]*Vovov[k,c,j,b]
        P_OoVv[i,j,a,b] -= T1[i,c]*T1[k,b]*Voovv[j,k,c,a]
        P_OoVv[i,j,a,b] += 2.0*T2[i,k,a,c]*Vovov[k,c,j,b]
        P_OoVv[i,j,a,b] -= T2[i,k,a,c]*Voovv[j,k,c,b]
        P_OoVv[i,j,a,b] -= T2[k,j,a,c]*Voovv[i,k,c,b]
        P_OoVv[i,j,a,b] += -2.0*T1[l,b]*T2[i,k,a,c]*Vooov[l,j,k,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[k,i,a,c]*Vooov[l,j,k,c]
        P_OoVv[i,j,a,b] -= T1[j,c]*T2[i,k,d,b]*Vovvv[k,c,a,d]
        P_OoVv[i,j,a,b] -= T1[j,c]*T2[k,i,a,d]*Vovvv[k,d,b,c]
        P_OoVv[i,j,a,b] -= T1[j,c]*T2[i,k,a,d]*Vovvv[k,c,b,d]
        P_OoVv[i,j,a,b] += T1[j,c]*T2[l,k,a,b]*Vooov[l,i,k,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[i,k,a,c]*Vooov[k,j,l,c]
        P_OoVv[i,j,a,b] -= T1[k,a]*T2[i,j,d,c]*Vovvv[k,d,b,c]
        P_OoVv[i,j,a,b] += T1[k,a]*T2[i,l,c,b]*Vooov[l,j,k,c]
        P_OoVv[i,j,a,b] += 2.0*T1[j,c]*T2[i,k,a,d]*Vovvv[k,d,b,c]
        P_OoVv[i,j,a,b] -= T1[k,c]*T2[i,j,a,d]*Vovvv[k,d,b,c]
        P_OoVv[i,j,a,b] += 2.0*T1[k,c]*T2[i,j,a,d]*Vovvv[k,c,b,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T2[i,l,a,b]*Vooov[k,j,l,c]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T2[i,l,a,b]*Vooov[l,j,k,c]
        P_OoVv[i,j,a,b] += T2[j,k,c,d]*T2[i,l,a,b]*Vovov[k,c,l,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[j,d]*T2[i,l,a,b]*Vovov[k,c,l,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[j,d]*T2[i,l,a,b]*Vovov[l,c,k,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[l,a]*T2[i,j,d,b]*Vovov[k,c,l,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[l,a]*T2[i,j,d,b]*Vovov[l,c,k,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,b,d]*Vovov[k,c,l,d]
        P_OoVv[i,j,a,b] += -2.0*T1[i,c]*T1[k,a]*T2[j,l,b,d]*Vovov[k,c,l,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,d,b]*Vovov[l,c,k,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[l,b]*T2[k,j,a,d]*Vovov[k,c,l,d]
        P_OoVv[i,j,a,b] += -2.0*T2[i,k,d,c]*T2[l,j,a,b]*Vovov[k,c,l,d]
     
        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end
end
