function mo_from_ao!(I::IntegralHelper{T,E,O}, aoints::IntegralHelper{T,E,AtomicOrbitals}, entries...) where {T<:AbstractFloat,E<:AbstractERI,
                                                                                                              O<:AbstractRestrictedOrbitals}
    t = @elapsed begin
        output("Converting integrals:")
        output(replace("($T, $E, AtomicOrbitals) ⇒ ($T, $E, $O)", "Fermi.Orbitals."=>""))
        for entry in entries
            output("• {:6s}", entry, ending=" ")
            compute!(I, entry, aoints)
            output("✔️")
        end
    end
    output("Done in {:5.5f} seconds.\n", t)
end

function mo_from_ao!(I::IntegralHelper{T1,E1,O}, aoints::IntegralHelper{T2,E2,AtomicOrbitals}, entries...) where {T1<:AbstractFloat,E1<:AbstractERI,T2<:AbstractFloat,
                                                                                                              E2<:AbstractERI,O<:AbstractRestrictedOrbitals}
    t = @elapsed begin
        output("Converting integrals:")
        output(replace("($T2, $E2, AtomicOrbitals) ⇒ ($T1, $E1, $O)", "Fermi.Orbitals."=>""))
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        newaoints = IntegralHelper{T1}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        output("IntegralHelper have different ERI handlers. A new AO object will be used:", ending=" ")
        output("{}", typeof(newaoints))

        for entry in entries
            if occursin(r"F[dijab]{0,2}", entry)
                output("Fock matrix will be computed using original ERI")
                output("• {:6s}", entry, ending=" ")
                compute_F!(I, aoints)
                output("✔️")
            end
        end
        entries = [entries...]
        filter!(i->!occursin(r"F[dijab]{0,2}",i), entries)

        for entry in entries
            output("• {:6s}", entry, ending=" ")
            compute!(I, entry, newaoints)
            output("✔️")
        end
    end
    output("Done in {:5.5f} seconds.\n", t)
end

# function mo_from_ao!(I::IntegralHelper{T1,E1,O}, aoints::IntegralHelper{T2,E2,AtomicOrbitals}, entries...) where {T1<:AbstractFloat,T2<:AbstractFloat,
#                                                             E1<:AbstractERI,E2<:AbstractERI,O<:AbstractRestrictedOrbitals}
#     if T1 !== T2 || ((E2 <: AbstractDFERI) ⊻ (E1 <: AbstractDFERI))
#         output("!! AO Integrals are not the same type as the MO. New integrals will be computed.")
#         for entry in entries
#             if occursin(r"F[dijab]{0,2}", entry)
#                 output("Fock matrix will be computed using old ERI")
#                 compute_F!(I, aoints)
#             end
#             entries = [entries...]
#             filter!(i->i==entry, entries)
#         end
#         basis = I.orbitals.basis
#         aoorbs = AtomicOrbitals(I.molecule, basis)
#         aoints = IntegralHelper{T1}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
#     end
#     t = @elapsed begin
#         output("Computing MO Integrals...")
#         for entry in entries
#             compute!(I, entry, aoints)
#         end
#     end
#     output("Done in {:5.5f} seconds.", t)
# end

# function mo_from_ao!(I::IntegralHelper{T1,E1,O}, aoints::IntegralHelper{T2,SparseERI,AtomicOrbitals}, entries...) where {T1<:AbstractFloat,T2<:AbstractFloat,
#                                                             E1<:AbstractERI,O<:AbstractRestrictedOrbitals}
#     if T1 !== T2 || E1 <: AbstractDFERI
#         output("!! AO Integrals are not the same type as the MO. New integrals will be computed.")
#         for entry in entries
#             if occursin(r"F[dijab]{0,2}", entry)
#                 output("Fock matrix will be computed using old ERI")
#                 compute_F!(I, aoints)
#             end
#             entries = [entries...]
#             filter!(i->i==entry, entries)
#         end
#         basis = I.orbitals.basis
#         aoorbs = AtomicOrbitals(I.molecule, basis)
#         aoints = IntegralHelper{T1}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
#     end
#     t = @elapsed begin
#         output("Converting integrals:")
#         output(replace("($T2, SparseERI, AtomicOrbitals) ⇒ ($T1, $E1, $O)", "Fermi.Orbitals."=>""))
#         for entry in entries
#             compute!(I, entry, aoints)
#         end
#     end
#     output("Done in {:5.5f} seconds.\n", t)
# end


function compute!(I::IntegralHelper{T,Chonky,O}, entry::String, x...) where {T<: AbstractFloat,O<:AbstractRestrictedOrbitals}
    if entry == "S"
        compute_S!(I, x...)
    elseif entry == "T"
        compute_T!(I, x...)
    elseif entry == "V"
        compute_V!(I, x...)
    elseif entry == "ERI"
        compute_ERI!(I, x...)
    elseif entry == "OOOO"
        compute_OOOO!(I, x...)
    elseif entry == "OOOV"
        compute_OOOV!(I, x...)
    elseif entry == "OVOV"
        compute_OVOV!(I, x...)
    elseif entry == "OOVV"
        compute_OOVV!(I, x...)
    elseif entry == "OVVV"
        compute_OVVV!(I, x...)
    elseif entry == "VVVV"
        compute_VVVV!(I, x...)
    elseif occursin(r"F[dijab]{0,2}", entry)
        compute_F!(I, x...)
    else
        throw(FermiException("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute!(I::IntegralHelper{T,E,O}, entry::String, x...) where {T<: AbstractFloat, E<:AbstractDFERI, O<:AbstractRestrictedOrbitals}
    if entry == "S"
        compute_S!(I, x...)
    elseif entry == "T"
        compute_T!(I, x...)
    elseif entry == "V"
        compute_V!(I, x...)
    elseif entry == "ERI"
        compute_ERI!(I, x...)
    elseif entry == "BOO"
        compute_BOO!(I, x...)
    elseif entry == "BOV"
        compute_BOV!(I, x...)
    elseif entry == "BVV"
        compute_BVV!(I, x...)
    elseif entry == "OOOO"
        compute_OOOO!(I, x...)
    elseif entry == "OOOV"
        compute_OOOV!(I, x...)
    elseif entry == "OVOV"
        compute_OVOV!(I, x...)
    elseif entry == "OOVV"
        compute_OOVV!(I, x...)
    elseif entry == "OVVV"
        compute_OVVV!(I, x...)
    elseif entry == "VVVV"
        compute_VVVV!(I, x...)
    elseif occursin(r"F[dijab]{0,2}", entry)
        compute_F!(I, x...)
    else
        throw(FermiException("Invalid key for IntegralHelper: $(entry)."))
    end
end

function compute_ERI!(I::IntegralHelper{T, E, O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
        # Create AO integral object
        basis = I.orbitals.basis
        aoorbs = AtomicOrbitals(I.molecule, basis)
        aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
        compute_ERI!(I,aoints)
end

include("ROIntegrals/OneElectron.jl")
include("ROIntegrals/DFERI.jl")

function compute_OOOO!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Create AO integral object
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_OOOO!(I,aoints)
end

function compute_OOOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Create AO integral object
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_OOOV!(I,aoints)
end

function compute_OOVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Create AO integral object
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_OOVV!(I,aoints)
end

function compute_OVOV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Create AO integral object
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_OVOV!(I,aoints)
end

function compute_OVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Create AO integral object
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_OVVV!(I,aoints)
end
function compute_VVVV!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    # Create AO integral object
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper{T}(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_VVVV!(I,aoints)
end


function compute_ERI!(I::IntegralHelper{T, Chonky, O}, aoints::IntegralHelper{T, Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
    AOERI = aoints["ERI"]
    C = I.orbitals.C
    if eltype(I.orbitals.C) !== T
        C = T.(C)
    end
    @tensoropt MOERI[i,j,k,l] :=  AOERI[μ, ν, ρ, σ]*C[μ, i]*C[ν, j]*C[ρ, k]*C[σ, l]
    I["ERI"] = MOERI
end

function compute_OOOO!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    core = Options.get("drop_occ")
    ndocc = I.molecule.Nα
    o = (1+core):ndocc
    Co = I.orbitals.C[:,o]
    if eltype(I.orbitals.C) !== T
        Co = T.(Co)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOOO[i,j,k,l] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Co[σ, l]
    end
    I["OOOO"] = OOOO
end

function compute_OOOO!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    eri_vals = aoints["ERI"].data
    idxs = aoints["ERI"].indexes

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
    @sync for i = o
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

function compute_OOOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

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
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOOV[i,j,k,a] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Co[ρ, k]*Cv[σ, a]
    end
    I["OOOV"] = OOOV
end

function compute_OOOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

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
    @sync for i = o
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

function compute_OOVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

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
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OOVV[i,j,a,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Co[ν, j]*Cv[ρ, a]*Cv[σ, b]
    end
    I["OOVV"] = OOVV
end

function compute_OOVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

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
    @sync for i = o
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

function compute_OVOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

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
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OVOV[i,a,j,b] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Co[ρ, j]*Cv[σ, b]
    end
    I["OVOV"] = OVOV
end

function compute_OVOV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

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
    @sync for i = o
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

function compute_OVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

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
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        OVVV[i,a,b,c] :=  AOERI[μ, ν, ρ, σ]*Co[μ, i]*Cv[ν, a]*Cv[ρ, b]*Cv[σ, c]
    end
    I["OVVV"] = OVVV
end

function compute_OVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

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
    @sync for i = o
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

function compute_VVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,Chonky, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

    AOERI = aoints["ERI"]

    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    v = (ndocc+1):(nbf - inac)
    Cv = I.orbitals.C[:,v]
    if eltype(I.orbitals.C) !== T
        Cv = T.(Cv)
    end
    @tensoropt (μ=>100x, ν=>100x, ρ=>100x, σ=>100x, i=>10x, j=>10x, k=>10x, l=>10x, a=>80x, b=>80x, c=>80x, d=>80) begin 
        VVVV[a,b,c,d] :=  AOERI[μ, ν, ρ, σ]*Cv[μ, a]*Cv[ν, b]*Cv[ρ, c]*Cv[σ, d]
    end
    I["VVVV"] = VVVV
end

function compute_VVVV!(I::IntegralHelper{T,Chonky,O}, aoints::IntegralHelper{T,SparseERI, AtomicOrbitals}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}

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

function compute_F!(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}
    basis = I.orbitals.basis
    aoorbs = AtomicOrbitals(I.molecule, basis)
    aoints = IntegralHelper(molecule=I.molecule, orbitals=aoorbs, basis=basis, eri_type=I.eri_type)
    compute_F!(I, aoints)
end

function compute_F!(I::IntegralHelper{<:AbstractFloat,Chonky,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{<:AbstractFloat, Chonky, AtomicOrbitals})

    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)

    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)

    Co = I.orbitals.C[:,1:ndocc]
    if eltype(I.orbitals.C) !== T
        Co = T.(Cv)
    end

    # Build Atomic Fock
    F = aoints["T"] + aoints["V"]
    ERI = aoints["ERI"]
    D = Co*transpose(Co)
    
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
    C = I.orbitals.C

    # Atomic Fock -> Molecular Fock
    @tensoropt Fmol[p,q] := C[μ,p]*F[μ,ν]*C[ν,q]

    Fd = FermiMDArray(diag(Fmol))
    I["Fd"] = Fd
    I["Fii"] = Fd[o]
    I["Faa"] = Fd[v]
    I["Fia"] = Fmol[o,v] 
    I["Fij"] = Fmol[o,o] - diagm(Fd[o])
    I["Fab"] = Fmol[v,v] - diagm(Fd[v])
end

function compute_F!(I::IntegralHelper{T1,E1,<:AbstractRestrictedOrbitals}, aoints::IntegralHelper{T2,E2,AtomicOrbitals}) where {T1<:AbstractFloat,T2<:AbstractFloat,
                                                                                                                              E1<:AbstractDFERI,E2<:AbstractDFERI}

    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)

    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)

    Co = I.orbitals.C[:,1:ndocc]
    if eltype(I.orbitals.C) !== T1
        Co = T1.(Co)
    end

    # Build Atomic Fock
    F = aoints["T"] + aoints["V"]
    b = aoints["ERI"]
    D = Co*transpose(Co)
    
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
    C = I.orbitals.C

    # Atomic Fock -> Molecular Fock
    @tensoropt Fmol[p,q] := C[μ,p]*F[μ,ν]*C[ν,q]

    Fd = FermiMDArray(diag(Fmol))
    I["Fd"] = Fd
    I["Fii"] = Fd[o]
    I["Faa"] = Fd[v]
    I["Fia"] = Fmol[o,v] 
    I["Fij"] = Fmol[o,o] - diagm(Fd[o])
    I["Fab"] = Fmol[v,v] - diagm(Fd[v])
end

function compute_F!(I::IntegralHelper{<:AbstractFloat,Chonky,RHFOrbitals}, aoints::IntegralHelper{<:AbstractFloat, SparseERI, AtomicOrbitals})
    compute_F!(I)
end

function compute_F!(I::IntegralHelper{<:AbstractFloat,RIFIT,RHFOrbitals}, aoints::IntegralHelper{<:AbstractFloat, JKFIT, AtomicOrbitals})
    compute_F!(I)
end

function compute_F!(I::IntegralHelper{T,E,RHFOrbitals}) where {T<:AbstractFloat, E<:AbstractERI}
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")
    ndocc = I.molecule.Nα
    nbf = size(I.orbitals.C,1)
    o = (1+core):ndocc
    v = (ndocc+1):(nbf - inac)
    No = length(o)
    Nv = length(v)

    I["Fd"] = T.(I.orbitals.eps)
    I["Fii"] = T.(I.orbitals.eps[o])
    I["Faa"] = T.(I.orbitals.eps[v])
    I["Fia"] = FermiMDzeros(T, No,Nv)
    I["Fij"] = FermiMDzeros(T, No,No)
    I["Fab"] = FermiMDzeros(T, Nv,Nv)
end

function compute_ref_energy(I::IntegralHelper{T,Chonky,O}) where {T<:AbstractFloat, O<:AbstractRestrictedOrbitals}
    ndocc = I.molecule.Nα
    o = 1:ndocc

    H = I["T"] + I["V"]

    OOOO = I["OOOO"]

    E0 = zero(T)
    for i in o
        E0 += 2*H[i,i]
        for j in o
            E0 += 2*OOOO[i,i,j,j] - OOOO[i,j,i,j]
        end
    end

    return E0
end

function compute_ref_energy(I::IntegralHelper{T,E,O}) where {T<:AbstractFloat,E<:AbstractDFERI,O<:AbstractRestrictedOrbitals}
    ndocc = I.molecule.Nα
    o = 1:ndocc

    H = I["T"] + I["V"]

    Boo = I["BOO"]

    E0 = zero(T)
    for i in o
        E0 += 2*H[i,i]
        for j in o
            E0 += 2*sum(Boo[:,i,i] .* Boo[:,j,j])
            E0 -= sum(Boo[:,i,j] .* Boo[:,i,j])
        end
    end

    return E0
end