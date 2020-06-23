function form_Dia(oovv,F::Union{Array{Float64,2},Array{Float32,2}})
    dt = eltype(F)
    nocc = size(oovv)[1]
    nvir = size(oovv)[4]
    Dia = zeros(nocc,nvir)
    for i in 1:nocc
        for a in 1:nvir
            aa = a+nocc
            Dia[i,a] = F[i,i] - F[aa,aa]
        end
    end
    return Dia
end
function form_Dijab(tiJaB, F::Union{Array{Float64,1},Array{Float32,1}})
    dt = eltype(tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    Dijab = zeros(dt, nocc, nocc, nvir, nvir)
    for i in rocc
        for j in rocc
            for a in rvir
                for b in rvir
                    aa = a + nocc
                    bb = b + nocc
                    Dijab[i, j, a, b] = F[i] + F[j] - F[aa] - F[bb]
                end
            end
        end
    end
    return Dijab
end
function form_Dijab(ijab,f::Union{Array{Float64,2},Array{Float32,2}})
    nocc = size(ijab)[1]
    nvir = size(ijab)[4]
    Dijab = zeros(size(ijab))
    for i in 1:nocc
        for j in 1:nocc
            for a in 1:nvir
                for b in 1:nvir
                    aa = a + nocc
                    bb = b + nocc
                    Dijab[i,j,a,b] = f[i,i] + f[j,j] - f[aa,aa] - f[bb,bb]
                end
            end
        end
    end
    return Dijab
end
function form_DiJaB(iJaB,fa::Array{Float64,2},fb::Array{Float64,2})
    nocca = size(iJaB)[1]
    noccb = size(iJaB)[2]
    nvira = size(iJaB)[3]
    nvirb = size(iJaB)[4]
    DiJaB = zeros(size(iJaB))
    for i in 1:nocca
        for j in 1:noccb
            for a in 1:nvira
                for b in 1:nvirb
                    aa = a + nocca
                    bb = b + noccb
                    DiJaB[i,j,a,b] = fa[i,i] + fb[j,j] - fa[aa,aa] - fb[bb,bb]
                end
            end
        end
    end
    return DiJaB
end
