function RMP2{T}(refWfn::Wfn,alg::Fermi.MollerPlesset.Direct) where T <: AbstractFloat
    dmp2 = 0.0
    nocc = refWfn.nalpha
    rocc = 1:1:refWfn.nalpha
    epsa = refWfn.epsa
    Ca = refWfn.Ca
    mints = refWfn.mints
    basis = refWfn.basis
    nao = basis.nbf()
    nvir = nao - nocc
    rvir = 1:nvir
    eps = refWfn.epsa
    dmp2 = 0.0
    for i in rocc
        Cav = Ca[:,nocc+1:nao]
        Cao = Ca[:,1:nocc]
        D = zeros(nvir,nocc,nvir)
        for j in rocc
            for a in rvir
                for b in rvir
                    aa = a + nocc
                    bb = b + nocc
                    D[a,j,b] = 1/(eps[i] + eps[j] - eps[aa] - eps[bb])
                end
            end
        end
        temp = Direct.direct_ao_contract(i,Ca,mints,basis)
        @tensoropt begin
            temp2[a,λ,σ] := Cav[ν,a]*temp[ν,λ,σ]
        end
        temp = nothing
        @tensoropt begin
            temp3[a,j,σ] := Cao[λ,j]*temp2[a,λ,σ]
        end
        temp2 = nothing
        @tensoropt begin
            temp4[a,j,b] := Cav[σ,b]*temp3[a,j,σ]
        end
        temp3 = nothing
        out = temp4 .* (2*temp4 - permutedims(temp4,[3,2,1])) .* D
        dmp2 += reduce(+,out)
    end
    RMP2{T}(dmp2,Fermi.MemTensor{T}(zeros(T,0,0,0,0)),Fermi.MemTensor{T}(zeros(T,0,0,0,0)))
end

