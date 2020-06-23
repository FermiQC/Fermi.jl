"""
    do_direct_rmp2

Algorithm
--
for each occ <ij|-->, compute partially transformed
integrals <ij|λσ>, store in memory

for each vir <--|ab> contract <ij|λσ> ⋅ C[a,λ] ⋅ C[b,σ]
and its transpose <--|ba>

compute contribution of <ij|ab> to δE(MP2)

   12 34      12 34     12 34        
= <ij|ab> ( 2<ij|ab> - <ij|ba> ) / Δ[ijab]
= <ij|ab> ( 2<ij|ab> - <ij|ba> ) / Δ[ijab]
------------------------------------------
   13 24      13 24     13 24        
= [ia|jb] ( 2[ia|jb] - [ib|ja] ) / Δ[ijab]
"""
function do_direct_rmp2(refWfn::Wfn)
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
    return dmp2

end

