function cc_update_energy(T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,E,O}, 
                          alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Bov = moints["BOV"]
    @tensoropt (k=>x, l=>x, c=>100x, d=>100x, Q=>200x)  begin
        X[l,c,k,d] := 2.0*T2[k,l,c,d]
        X[l,c,k,d] -= T1[l,c]*T1[k,d]
        X[l,c,k,d] -= T2[l,k,c,d]
        CC_energy = X[l,c,k,d]*Bov[Q,k,c]*Bov[Q,l,d]
        CC_energy += 2.0*T1[l,c]*T1[k,d]*Bov[Q,l,c]*Bov[Q,k,d]
    end
    
    return CC_energy
end

function cc_update_T1!(newT1::AbstractArray{T,2}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,E,O}, 
                      alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Boo, Bov, Bvv = moints["BOO"], moints["BOV"], moints["BVV"]

    @tensoropt (i=>1, j=>1, k=>1, l=>1, a=>10, b=>10, c=>10, d=>10, Q=>20) begin
        newT1[i,a] -= T1[k,c]*Boo[Q,i,k]*Bvv[Q,c,a]
        newT1[i,a] += 2.0*T1[k,c]*Bov[Q,k,c]*Bov[Q,i,a]
        newT1[i,a] -= T2[k,i,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Boo[Q,k,i]*Bov[Q,l,c]
        newT1[i,a] += T2[l,k,a,c]*Boo[Q,k,i]*Bov[Q,l,c]
        newT1[i,a] += -2.0*T1[k,c]*T1[l,a]*Boo[Q,l,i]*Bov[Q,k,c]
        newT1[i,a] -= T1[k,c]*T1[i,d]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Bov[Q,k,c]*Bvv[Q,a,d]
        newT1[i,a] += T1[k,c]*T1[l,a]*Boo[Q,k,i]*Bov[Q,l,c]
        newT1[i,a] += -2.0*T1[k,c]*T2[i,l,a,d]*Bov[Q,l,c]*Bov[Q,k,d]
        newT1[i,a] += -2.0*T1[k,c]*T2[l,i,a,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT1[i,a] += T1[k,c]*T2[l,i,a,d]*Bov[Q,l,c]*Bov[Q,k,d]
        newT1[i,a] += -2.0*T1[i,c]*T2[l,k,a,d]*Bov[Q,l,c]*Bov[Q,k,d]
        newT1[i,a] += T1[i,c]*T2[l,k,a,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT1[i,a] += -2.0*T1[l,a]*T2[i,k,d,c]*Bov[Q,k,c]*Bov[Q,l,d]
        newT1[i,a] += T1[l,a]*T2[i,k,c,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT1[i,a] += T1[k,c]*T1[i,d]*T1[l,a]*Bov[Q,l,c]*Bov[Q,k,d]
        newT1[i,a] += -2.0*T1[k,c]*T1[i,d]*T1[l,a]*Bov[Q,k,c]*Bov[Q,l,d]
        newT1[i,a] += 4.0*T1[k,c]*T2[i,l,a,d]*Bov[Q,k,c]*Bov[Q,l,d]
    end
end

function cc_update_T2!(newT2::AbstractArray{T,4}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,E,O}, 
                    alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}

    Boo, Bov, Bvv = moints["BOO"], moints["BOV"], moints["BVV"]

    @tensoropt (i=>1, j=>1, k=>1, l=>1, a=>10, b=>10, c=>10, d=>10, Q=>20) begin
        newT2[i,j,a,b] += Bov[Q,i,a]*Bov[Q,j,b]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*Bvv[Q,c,a]*Bvv[Q,d,b]
        newT2[i,j,a,b] += T2[i,j,c,d]*Bvv[Q,c,a]*Bvv[Q,d,b]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*Boo[Q,i,k]*Boo[Q,j,l]
        newT2[i,j,a,b] += T2[k,l,a,b]*Boo[Q,i,k]*Boo[Q,j,l]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,a]*Bov[Q,k,c]*Bvv[Q,b,d]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,b]*Bov[Q,k,d]*Bvv[Q,a,c]
        newT2[i,j,a,b] += T1[i,c]*T1[k,a]*T1[l,b]*Boo[Q,l,j]*Bov[Q,k,c]
        newT2[i,j,a,b] += T1[j,c]*T1[k,a]*T1[l,b]*Boo[Q,k,i]*Bov[Q,l,c]
        newT2[i,j,a,b] += T2[k,l,a,c]*T2[i,j,d,b]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[l,j,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += -2.0*T2[l,k,a,c]*T2[i,j,d,b]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,d,b]*Bov[Q,l,c]*Bov[Q,k,d]
        newT2[i,j,a,b] += T2[i,k,a,c]*T2[l,j,b,d]*Bov[Q,l,c]*Bov[Q,k,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[j,l,b,d]*Bov[Q,l,c]*Bov[Q,k,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += -2.0*T2[k,i,a,c]*T2[j,l,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += T2[i,j,a,c]*T2[l,k,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += -2.0*T2[i,j,a,c]*T2[k,l,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += T2[k,j,a,c]*T2[i,l,d,b]*Bov[Q,l,c]*Bov[Q,k,d]
        newT2[i,j,a,b] += 4.0*T2[i,k,a,c]*T2[j,l,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += T2[i,j,d,c]*T2[l,k,a,b]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T1[k,a]*T1[l,b]*Bov[Q,k,c]*Bov[Q,l,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T2[l,k,a,b]*Bov[Q,l,c]*Bov[Q,k,d]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*T2[i,j,d,c]*Bov[Q,l,c]*Bov[Q,k,d]
        P_OoVv[i,j,a,b] := T1[j,c]*Bov[Q,i,a]*Bvv[Q,c,b]
        P_OoVv[i,j,a,b] -= T1[k,b]*Boo[Q,j,k]*Bov[Q,i,a]
        P_OoVv[i,j,a,b] -= T2[k,i,a,c]*Bov[Q,k,c]*Bov[Q,j,b]
        P_OoVv[i,j,a,b] -= T1[i,c]*T1[k,a]*Bov[Q,k,c]*Bov[Q,j,b]
        P_OoVv[i,j,a,b] -= T1[i,c]*T1[k,b]*Boo[Q,j,k]*Bvv[Q,c,a]
        P_OoVv[i,j,a,b] += 2.0*T2[i,k,a,c]*Bov[Q,k,c]*Bov[Q,j,b]
        P_OoVv[i,j,a,b] -= T2[i,k,a,c]*Boo[Q,j,k]*Bvv[Q,c,b]
        P_OoVv[i,j,a,b] -= T2[k,j,a,c]*Boo[Q,i,k]*Bvv[Q,c,b]
        P_OoVv[i,j,a,b] += -2.0*T1[l,b]*T2[i,k,a,c]*Boo[Q,l,j]*Bov[Q,k,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[k,i,a,c]*Boo[Q,l,j]*Bov[Q,k,c]
        P_OoVv[i,j,a,b] -= T1[j,c]*T2[i,k,d,b]*Bov[Q,k,c]*Bvv[Q,a,d]
        P_OoVv[i,j,a,b] -= T1[j,c]*T2[k,i,a,d]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] -= T1[j,c]*T2[i,k,a,d]*Bov[Q,k,c]*Bvv[Q,b,d]
        P_OoVv[i,j,a,b] += T1[j,c]*T2[l,k,a,b]*Boo[Q,l,i]*Bov[Q,k,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[i,k,a,c]*Boo[Q,k,j]*Bov[Q,l,c]
        P_OoVv[i,j,a,b] -= T1[k,a]*T2[i,j,d,c]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += T1[k,a]*T2[i,l,c,b]*Boo[Q,l,j]*Bov[Q,k,c]
        P_OoVv[i,j,a,b] += 2.0*T1[j,c]*T2[i,k,a,d]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] -= T1[k,c]*T2[i,j,a,d]*Bov[Q,k,d]*Bvv[Q,b,c]
        P_OoVv[i,j,a,b] += 2.0*T1[k,c]*T2[i,j,a,d]*Bov[Q,k,c]*Bvv[Q,b,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T2[i,l,a,b]*Boo[Q,k,j]*Bov[Q,l,c]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T2[i,l,a,b]*Boo[Q,l,j]*Bov[Q,k,c]
        P_OoVv[i,j,a,b] += T2[j,k,c,d]*T2[i,l,a,b]*Bov[Q,k,c]*Bov[Q,l,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[j,d]*T2[i,l,a,b]*Bov[Q,k,c]*Bov[Q,l,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[j,d]*T2[i,l,a,b]*Bov[Q,l,c]*Bov[Q,k,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[l,a]*T2[i,j,d,b]*Bov[Q,k,c]*Bov[Q,l,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[l,a]*T2[i,j,d,b]*Bov[Q,l,c]*Bov[Q,k,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        P_OoVv[i,j,a,b] += -2.0*T1[i,c]*T1[k,a]*T2[j,l,b,d]*Bov[Q,k,c]*Bov[Q,l,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,d,b]*Bov[Q,l,c]*Bov[Q,k,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[l,b]*T2[k,j,a,d]*Bov[Q,k,c]*Bov[Q,l,d]
        P_OoVv[i,j,a,b] += -2.0*T2[i,k,d,c]*T2[l,j,a,b]*Bov[Q,k,c]*Bov[Q,l,d]

        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end
end
