function compute_OOOO!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    idxs = aoints["ERI"].indexes
    eri_vals = aoints["ERI"].data

    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Co = T.(Co)
    end

    # Create a partially-contracted dense array from SparseERI
    iνρσ = zeros(T, length(o), nbf, nbf, nbf)
    @sync for i = 1:length(o)
        Threads.@spawn begin
        @inbounds begin
        @fastmath begin
        for z = eachindex(eri_vals)
            V = eri_vals[z]
            μ,ν,ρ,σ = idxs[z] .+ 1
            γμν = μ !== ν
            γρσ = ρ !== σ
            γab = Fermi.index2(μ,ν) !== Fermi.index2(ρ, σ)
            Vμ = V*Co[μ,i] 
            Vν = V*Co[ν,i] 
            Vρ = V*Co[ρ,i] 
            Vσ = V*Co[σ,i] 

            if γab && γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
                iνρσ[i,ρ,ν,μ] += Vσ 

            elseif γab && γμν
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 

            elseif γab && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
            
            elseif γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
            elseif γab
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
            else
                iνρσ[i,ν,ρ,σ] += Vμ 
            end
        end
        end
        end
        end
    end

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOOO[i,j,k,l] :=  iνρσ[i, ν, ρ, σ]*Co[ν, j]*Co[ρ, k]*Co[σ, l]
    end
    I["OOOO"] = OOOO
end


function compute_OOOV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]

    # Create a partially-contracted dense array from SparseERI
    iνρσ = zeros(T, length(o), nbf, nbf, nbf)
    @sync for i = 1:length(o)
        Threads.@spawn begin
        @inbounds begin
        @fastmath begin
        for z = eachindex(eri_vals)
            V = eri_vals[z]
            μ,ν,ρ,σ = idxs[z] .+ 1
            γμν = μ !== ν
            γρσ = ρ !== σ
            γab = Fermi.index2(μ,ν) !== Fermi.index2(ρ, σ)
            Vμ = V*Co[μ,i] 
            Vν = V*Co[ν,i] 
            Vρ = V*Co[ρ,i] 
            Vσ = V*Co[σ,i] 

            if γab && γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
                iνρσ[i,ρ,ν,μ] += Vσ 

            elseif γab && γμν
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 

            elseif γab && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
            
            elseif γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
            elseif γab
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
            else
                iνρσ[i,ν,ρ,σ] += Vμ 
            end
        end
        end
        end
        end
    end

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOOV[i,j,k,a] :=  iνρσ[i, ν, ρ, σ]*Co[ν, j]*Co[ρ, k]*Cv[σ, a]
    end
    I["OOOV"] = OOOV
end


function compute_OOVV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end

    # Create a partially-contracted dense array from SparseERI
    iνρσ = zeros(T, length(o), nbf, nbf, nbf)
    @sync for i = 1:length(o)
        Threads.@spawn begin
        @inbounds begin
        @fastmath begin
        for z = eachindex(eri_vals)
            V = eri_vals[z]
            μ,ν,ρ,σ = idxs[z] .+ 1
            γμν = μ !== ν
            γρσ = ρ !== σ
            γab = Fermi.index2(μ,ν) !== Fermi.index2(ρ, σ)
            Vμ = V*Co[μ,i] 
            Vν = V*Co[ν,i] 
            Vρ = V*Co[ρ,i] 
            Vσ = V*Co[σ,i] 

            if γab && γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
                iνρσ[i,ρ,ν,μ] += Vσ 

            elseif γab && γμν
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 

            elseif γab && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
            
            elseif γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
            elseif γab
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
            else
                iνρσ[i,ν,ρ,σ] += Vμ 
            end
        end
        end
        end
        end
    end

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOVV[i,j,a,b] :=  iνρσ[i, ν, ρ, σ]*Co[ν, j]*Cv[ρ, a]*Cv[σ, b]
    end
    I["OOVV"] = OOVV
end


function compute_OVOV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end

    # Create a partially-contracted dense array from SparseERI
    iνρσ = zeros(T, length(o), nbf, nbf, nbf)
    @sync for i = 1:length(o)
        Threads.@spawn begin
        @inbounds begin
        @fastmath begin
        for z = eachindex(eri_vals)
            V = eri_vals[z]
            μ,ν,ρ,σ = idxs[z] .+ 1
            γμν = μ !== ν
            γρσ = ρ !== σ
            γab = Fermi.index2(μ,ν) !== Fermi.index2(ρ, σ)
            Vμ = V*Co[μ,i] 
            Vν = V*Co[ν,i] 
            Vρ = V*Co[ρ,i] 
            Vσ = V*Co[σ,i] 

            if γab && γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
                iνρσ[i,ρ,ν,μ] += Vσ 

            elseif γab && γμν
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 

            elseif γab && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
            
            elseif γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
            elseif γab
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
            else
                iνρσ[i,ν,ρ,σ] += Vμ 
            end
        end
    end
end
    end
    end

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OVOV[i,a,j,b] :=  iνρσ[i, ν, ρ, σ]*Cv[ν, a]*Co[ρ, j]*Cv[σ, b]
    end
    I["OVOV"] = OVOV
end

function compute_OVVV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
        Co = T.(Co)
    end

    # Create a partially-contracted dense array from SparseERI
    iνρσ = zeros(T, length(o), nbf, nbf, nbf)
    @sync for i = 1:length(o)
        Threads.@spawn begin
        @inbounds begin
        @fastmath begin
        for z = eachindex(eri_vals)
            V = eri_vals[z]
            μ,ν,ρ,σ = idxs[z] .+ 1
            γμν = μ !== ν
            γρσ = ρ !== σ
            γab = Fermi.index2(μ,ν) !== Fermi.index2(ρ, σ)
            Vμ = V*Co[μ,i] 
            Vν = V*Co[ν,i] 
            Vρ = V*Co[ρ,i] 
            Vσ = V*Co[σ,i] 

            if γab && γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
                iνρσ[i,ρ,ν,μ] += Vσ 

            elseif γab && γμν
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,σ,ν,μ] += Vρ 

            elseif γab && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
                iνρσ[i,ρ,μ,ν] += Vσ 
            
            elseif γμν && γρσ
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,ν,σ,ρ] += Vμ 
                iνρσ[i,μ,ρ,σ] += Vν 
                iνρσ[i,μ,σ,ρ] += Vν 
            elseif γab
                iνρσ[i,ν,ρ,σ] += Vμ 
                iνρσ[i,σ,μ,ν] += Vρ 
            else
                iνρσ[i,ν,ρ,σ] += Vμ 
            end
        end
        end
        end
        end
    end

    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OVVV[i,a,b,c] := iνρσ[i, ν, ρ, σ]*Cv[ν, a]*Cv[ρ, b]*Cv[σ, c]
    end
    I["OVVV"] = OVVV
end


function compute_VVVV!(I::IntegralHelper{T,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where T<:AbstractFloat

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
    end

    # Create a partially-contracted dense array from SparseERI
    aνρσ = zeros(T, length(v), nbf, nbf, nbf)
    @sync for a = 1:length(v)
        Threads.@spawn begin
        @inbounds begin
        @fastmath begin
        for z = eachindex(eri_vals)
            V = eri_vals[z]
            μ,ν,ρ,σ = idxs[z] .+ 1
            γμν = μ !== ν
            γρσ = ρ !== σ
            γab = Fermi.index2(μ,ν) !== Fermi.index2(ρ, σ)
            Vμ = V*Cv[μ,a] 
            Vν = V*Cv[ν,a] 
            Vρ = V*Cv[ρ,a] 
            Vσ = V*Cv[σ,a] 

            if γab && γμν && γρσ
                aνρσ[a,ν,ρ,σ] += Vμ 
                aνρσ[a,ν,σ,ρ] += Vμ 
                aνρσ[a,μ,ρ,σ] += Vν 
                aνρσ[a,μ,σ,ρ] += Vν 
                aνρσ[a,σ,μ,ν] += Vρ 
                aνρσ[a,σ,ν,μ] += Vρ 
                aνρσ[a,ρ,μ,ν] += Vσ 
                aνρσ[a,ρ,ν,μ] += Vσ 

            elseif γab && γμν
                aνρσ[a,ν,ρ,σ] += Vμ 
                aνρσ[a,μ,ρ,σ] += Vν 
                aνρσ[a,σ,μ,ν] += Vρ 
                aνρσ[a,σ,ν,μ] += Vρ 

            elseif γab && γρσ
                aνρσ[a,ν,ρ,σ] += Vμ 
                aνρσ[a,ν,σ,ρ] += Vμ 
                aνρσ[a,σ,μ,ν] += Vρ 
                aνρσ[a,ρ,μ,ν] += Vσ 
            
            elseif γμν && γρσ
                aνρσ[a,ν,ρ,σ] += Vμ 
                aνρσ[a,ν,σ,ρ] += Vμ 
                aνρσ[a,μ,ρ,σ] += Vν 
                aνρσ[a,μ,σ,ρ] += Vν 
            elseif γab
                aνρσ[a,ν,ρ,σ] += Vμ 
                aνρσ[a,σ,μ,ν] += Vρ 
            else
                aνρσ[a,ν,ρ,σ] += Vμ 
            end
        end
        end
        end
        end
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        VVVV[a,b,c,d] :=  aνρσ[a, ν, ρ, σ]*Cv[ν, b]*Cv[ρ, c]*Cv[σ, d]
    end
    I["VVVV"] = VVVV
end