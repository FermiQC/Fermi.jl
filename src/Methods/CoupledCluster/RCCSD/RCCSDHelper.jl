"""
    Fermi.CoupledCluster.cc_update_energy(T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,O}, 
                          alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Compute CC energy from T1 and T2 amplitudes constructed with a set of restricted orbitals.
NOTE: It does not include off-diagonal Fock contributions: See `od_cc_update_energy`.
"""
function cc_update_energy(T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,Chonky,O}, 
                          alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    Voovv = moints["OOVV"]
    @tensoropt (k=>x, l=>x, c=>100x, d=>100x)  begin
        B[l,c,k,d] := 2.0*T2[k,l,c,d]
        B[l,c,k,d] -= T1[l,c]*T1[k,d]
        B[l,c,k,d] -= T2[l,k,c,d]
        CC_energy = B[l,c,k,d]*Voovv[k,l,c,d]
        CC_energy += 2.0*T1[l,c]*T1[k,d]*Voovv[l,k,c,d]
    end
    
    return CC_energy
end

"""
    Fermi.CoupledCluster.od_cc_update_energy(T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,O}, 
                          alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Compute contributions from off-diagonal Fock terms to the CC energy from T1 and T2 amplitudes constructed with a set of restricted orbitals.
See also: `cc_update_energy`
"""
function od_cc_update_energy(T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,Chonky,O}, 
                          alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    Voovv = moints["OOVV"]
    CC_energy = cc_update_energy(T1, T2, moints, alg)
    @tensor CC_energy += 2.0*f[k,c]*T1[k,c]
    
    return CC_energy
end
"""
    Fermi.CoupledCluster.cc_update_T1!(newT1::AbstractArray{T,2}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,O}, 
                      alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Compute DᵢₐTᵢₐ from old T1 and T2 amplitudes. Final updated T1 amplitudes can be obtained by applying the denominator 1/Dᵢₐ.
NOTE: It does not include off-diagonal Fock contributions: See `od_cc_update_T1`.
"""
function cc_update_T1!(newT1::AbstractArray{T,2}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,Chonky,O}, 
                      alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = moints["OOOO"], moints["OOOV"], moints["OOVV"], moints["OVOV"], moints["OVVV"], moints["VVVV"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
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

"""
    Fermi.CoupledCluster.od_cc_update_T1!(newT1::AbstractArray{T,2}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,O}, 
                      alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Compute non-diagonal Fock contributions to DᵢₐTᵢₐ from old T1 and T2 amplitudes. Final updated T1 amplitudes 
can be obtained by applying the denominator 1/Dᵢₐ.
See also: `cc_update_T1`.
"""
function od_cc_update_T1!(newT1::AbstractArray{T,2}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,E,O}, 
                      alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    # Include non-RHF terms
    fov = moints["Fia"]
    foo = moints["Fij"]
    fvv = moints["Fab"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
        newT1[i,a] += fov[i,a]
        newT1[i,a] -= foo[i,k]*T1[k,a]
        newT1[i,a] += fvv[c,a]*T1[i,c]
        newT1[i,a] -= fov[k,c]*T1[i,c]*T1[k,a]
        newT1[i,a] += 2.0*fov[k,c]*T2[i,k,a,c]
        newT1[i,a] -= fov[k,c]*T2[k,i,a,c]
    end

    # Include RHF terms
    cc_update_T1!(newT1, T1, T2, moints, alg)
end

"""
    Fermi.CoupledCluster.cc_update_T2!(newT2::AbstractArray{T,4}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,O}, 
                    alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Compute Dᵢⱼₐᵦ⋅Tᵢⱼₐᵦ from old T1 and T2 amplitudes. Final updated T2 amplitudes can be obtained by applying the denominator 1/Dᵢⱼₐᵦ.
NOTE: It does not include off-diagonal Fock contributions: See `od_cc_update_T2`.
"""
function cc_update_T2!(newT2::AbstractArray{T,4}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,Chonky,O}, 
                    alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = moints["OOOO"], moints["OOOV"], moints["OOVV"], moints["OVOV"], moints["OVVV"], moints["VVVV"]

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
        P_OoVv[i,j,a,b] := -1.0*T1[k,b]*Vooov[j,i,k,a]
        P_OoVv[i,j,a,b] += T1[j,c]*Vovvv[i,c,a,b]
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

"""
    Fermi.CoupledCluster.od_cc_update_T2!(newT2::AbstractArray{T,4}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,O}, 
                    alg::RCCSDa) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Compute non-diagonal Fock contribution to Dᵢⱼₐᵦ⋅Tᵢⱼₐᵦ from old T1 and T2 amplitudes. Final updated T2 amplitudes 
can be obtained by applying the denominator 1/Dᵢⱼₐᵦ. See also: `cc_update_T2`.
"""
function od_cc_update_T2!(newT2::AbstractArray{T,4}, T1::AbstractArray{T,2}, T2::AbstractArray{T,4}, moints::IntegralHelper{T,E,O}, 
                    alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    # Include non-RHF terms
    fov = moints["Fia"]
    foo = moints["Fij"]
    fvv = moints["Fab"]
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>10x, b=>10x, c=>10x, d=>10x) begin
        P_OoVv[i,j,a,b] := -1.0*foo[i,k]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += fvv[c,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] -= fov[k,c]*T1[i,c]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] -= fov[k,c]*T1[k,a]*T2[i,j,c,b]
        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end
    #Include RHF terms
    cc_update_T2!(newT2, T1, T2, moints, alg)
end

"""
    Fermi.CoupledCluster.update_amp!(newT1::AbstractArray{T,2}, newT2::Array{T,4}, T1::Array{T, 2}, T2::Array{T, 4}, moints::IntegralHelper{T,O}, 
                     alg::A) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Computes new T1 and T2 amplitudes from old ones. It assumes arbitrary restricted orbitals.
"""
function update_amp!(newT1::AbstractArray{T,2}, newT2::AbstractArray{T,4}, T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, 
                    moints::IntegralHelper{T,E,O}, alg::RCCSDa) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    # Clean the arrays
    fill!(newT1, 0.0)
    fill!(newT2, 0.0)

    # Get new amplitudes
    od_cc_update_T1!(newT1, T1, T2, moints, alg)
    od_cc_update_T2!(newT2, T1, T2, moints, alg)

    # Orbital energies line
    if haskey(moints.cache, "D1")
        d = moints["D1"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        ϵo = Fd[(1+Options.get("drop_occ"):ndocc)]
        ϵv = Fd[(1+ndocc):size(T1,2)]

        d = FermiMDArray([ϵo[i]-ϵv[a] for i=eachindex(ϵo), a=eachindex(ϵv)])
        moints["D1"] = d
    end

    if haskey(moints.cache, "D2")
        D = moints["D2"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        ϵo = Fd[(1+Options.get("drop_occ"):ndocc)]
        ϵv = Fd[(1+ndocc):size(T1,2)]

        D = FermiMDArray([ϵo[i]+ϵ[j]-ϵv[a]-ϵv[b] for i=eachindex(ϵo), j=eachindex(ϵo), a=eachindex(ϵv), b=eachindex(ϵv)])
        moints["D2"] = D
    end

    newT1 ./= d
    newT2 ./= D
end

"""
    Fermi.CoupledCluster.update_amp!(newT1::AbstractArray{T,2}, newT2::Array{T,4}, T1::Array{T, 2}, T2::Array{T, 4}, moints::IntegralHelper{T,O}, 
                     alg::A) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

Computes new T1 and T2 amplitudes from old ones. It assumes Restricted Hartree-Fock reference.
"""
function update_amp!(newT1::AbstractArray{T,2}, newT2::AbstractArray{T,4}, T1::AbstractArray{T, 2}, T2::AbstractArray{T, 4}, moints::IntegralHelper{T,E,RHFOrbitals}, 
                     alg::RCCSDa) where {T<:AbstractFloat,E<:AbstractERI}

    # Clean the arrays
    fill!(newT1, 0.0)
    fill!(newT2, 0.0)

    # Get new amplitudes
    cc_update_T1!(newT1, T1, T2, moints, alg)
    cc_update_T2!(newT2, T1, T2, moints, alg)

    # Orbital energies line
    if haskey(moints.cache, "D1")
        d = moints["D1"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        frozen = Options.get("drop_occ")
        inac = Options.get("drop_vir")
        ϵo = Fd[(1+frozen:ndocc)]
        ϵv = Fd[(1+ndocc):end-inac]

        d = FermiMDArray([ϵo[i]-ϵv[a] for i=eachindex(ϵo), a=eachindex(ϵv)])
        moints["D1"] = d
    end

    if haskey(moints.cache, "D2")
        D = moints["D2"]
    else
        Fd = moints["Fd"]
        ndocc = moints.molecule.Nα
        frozen = Options.get("drop_occ")
        inac = Options.get("drop_vir")
        ϵo = Fd[(1+frozen:ndocc)]
        ϵv = Fd[(1+ndocc):end-inac]

        D = FermiMDArray([ϵo[i]+ϵo[j]-ϵv[a]-ϵv[b] for i=eachindex(ϵo), j=eachindex(ϵo), a=eachindex(ϵv), b=eachindex(ϵv)])
        moints["D2"] = D
    end

    # Orbital energies line
    d, D = moints["D1"], moints["D2"]

    newT1 ./= d
    newT2 ./= D
end