using Fermi
using Fermi.Integrals
using TensorOperations
using Combinatorics
using SparseArrays
using ArnoldiMethod

function RFCI(alg::SparseHamiltonian)
    aoints = IntegralHelper{Float64}()
    rhf = Fermi.HartreeFock.RHF(aoints)
    RFCI(aoints, rhf, alg)
end

function RFCI(aoints::IntegralHelper, rhf::Fermi.HartreeFock.RHF, alg::SparseHamiltonian)
    ci_header()

    if typeof(aoints.eri_type) === JKFIT
        aoints = IntegralHelper(eri_type=RIFIT())
    elseif Options.get("precision") == "single"
        aoints = IntegralHelper()
    end
    moints = IntegralHelper(orbitals=rhf.orbitals)

    Eref = moints.orbitals.sd_energy
    mol = moints.molecule
    Vnuc = Molecules.nuclear_repulsion(mol.atoms)
    Nelec = mol.Nα
    Nbas = aoints.orbitals.basisset.nbas
    Nvir = Nbas - Nelec

    Nfrozen = Options.get("drop_occ")
    Ninac = Options.get("drop_vir")
    Nactive = Nbas - Nfrozen - Ninac

    act_elec = 2*(Nactive - Nvir + Ninac)

    if act_elec < 0
        error("\nInvalid number of frozen orbitals ($Nfrozen) for $(2*Nelec) electrons.")
    end

    if Nactive ≤ act_elec/2
        error("\nNumber of active orbitals ($Nactive) too small for $(act_elec) active electrons")
    end

    if Nactive+Ninac+Nfrozen > Nbas
        error("\nSum of active ($Nactive) frozen ($Nfrozen) and inactive ($Ninac) orbitals greater than number of orbitals ($Nbas)")
    end

    # Get integrals
    output("Transforming integrals... ", ending="")
    t = @elapsed begin
        hp  = Fermi.Integrals.compute!(moints, aoints, "T")
        hp += Fermi.Integrals.compute!(moints, aoints, "V")
        eri = Fermi.Integrals.compute!(moints, "ERI")
    end

    r = 1:(Nfrozen+Nactive)
    hp = hp[r,r]
    eri = eri[r,r,r,r]
    aoints = nothing
    moints = nothing
    output("Done in {:10.5f} seconds.\n", t)

    output(" => Active Space Information ({:d}e, {:d}o)", act_elec, Nactive)
    output(" • # of Total Electrons:       {:>5d}", Nelec)
    output(" • # of Active Electrons:      {:>5d}", act_elec)
    output(" • # of Orbitals:              {:>5d}", Nbas)
    output(" • # of Frozen Orbitals:       {:>5d}", Nfrozen)
    output(" • # of Inactive Orbitals:     {:>5d}", Ninac)
    output(" • # of Active Orbitals:       {:>5d}", Nactive)

    output("\nGenerating determinants ...")
    t = @elapsed begin
        dets = get_determinants(act_elec, Nactive, Nfrozen)
        Ndets = length(dets)
    end
    output(" done in {} seconds.\n", t)
    output("\nNumber of Determinants: {:10d}\n", Ndets)

    output("\nBuilding Sparse Hamiltonian...")
    t = @elapsed begin
        H = get_sparse_hamiltonian_matrix(dets, hp, eri, Fermi.Options.get("cas_cutoff"))
    end
    output(" done in {:5.5f} seconds.\n", t)

    nroot = Fermi.Options.get("cas_nroot")

    output("Hamiltonian Matrix size: {:10.3f} Mb\n", Base.summarysize(H)/10^6)
    output("Diagonalizing Hamiltonian for {:3d} eigenvalues...", nroot)
    t = @elapsed begin
        decomp, history = partialschur(H, nev=nroot, tol=10^-12, which=LM())
        λ, ϕ = partialeigen(decomp)
    end
    output(" done in {:5.5f} seconds.\n", t)
    Etotal = λ[1]+Vnuc
    output("\n Final FCI Energy: {:15.10f}\n", Etotal)

    # Sort dets by importance
    C = ϕ[:,1]
    sp = sortperm(abs.(ϕ[:,1]), rev=true)
    C = ϕ[:,1][sp]
    dets = dets[sp]

    output("\n • Most important configurations:\n\n")
    output("Coefficient / Determinant / α-Occupancy / β-Occupancy")
    for i in 1:(min(10,length(C)))
        output("{:15.5f}      {}", C[i], detstring(dets[i], Nfrozen+Nactive))
    end
    output("\n")

    return RFCI(Etotal, Etotal-Eref, C, dets)
end

function get_sparse_hamiltonian_matrix(dets::Vector{Determinant{D}}, h::Array{T,2}, V::Array{T,4}, tol::Float64) where {T <: AbstractFloat, D <: Integer}

    Ndets = length(dets)
    Nα = sum(αlist(dets[1]))
    Nβ = sum(αlist(dets[1]))

    αind = [Array{Int64,1}(undef,Nα) for i = 1:Threads.nthreads()]
    βind = [Array{Int64,1}(undef,Nβ) for i = 1:Threads.nthreads()]
    vals = [T[] for i = 1:Threads.nthreads()]
    ivals = [Int64[] for i = 1:Threads.nthreads()]
    jvals = [Int64[] for i = 1:Threads.nthreads()]

    @sync for i in 1:Ndets
        Threads.@spawn begin
        D1 = dets[i]
        αindex!(D1, αind[Threads.threadid()])
        βindex!(D1, βind[Threads.threadid()])
        elem = 0.0
        for j in i:Ndets
            D2 = dets[j]
            αexc = αexcitation_level(D1,D2)
            βexc = βexcitation_level(D1,D2)
            el = αexc + βexc
            if el > 2
                continue 
            elseif el == 2
                elem = Hd2(D1, D2, V, αexc)
            elseif el == 1
                elem = Hd1(αind[Threads.threadid()], βind[Threads.threadid()], D1, D2, h, V, αexc)
            else
                elem = Hd0(αind[Threads.threadid()], βind[Threads.threadid()], h, V)
            end
            if abs(elem) > tol
                push!(vals[Threads.threadid()], elem)
                push!(ivals[Threads.threadid()], i)
                push!(jvals[Threads.threadid()], j)
            end
        end
    end #Threads.@spawn
    end

    ivals = vcat(ivals...)
    jvals = vcat(jvals...)
    vals  = vcat(vals...)
    return Symmetric(sparse(ivals, jvals, vals))
end

function get_1p_coupling_coefficients(dets::Array{Determinant,1}, nmo::Int)

    Ndets = length(dets)
    Nα = sum(αlist(dets[1]))
    Nβ = sum(αlist(dets[1]))
    αind = Array{Int64,1}(undef,Nα)
    βind = Array{Int64,1}(undef,Nβ)    
    γ = zeros(nmo, nmo, Ndets, Ndets)

    for nI in 1:Ndets
        I = dets[nI]
        for nJ in 1:Ndets
            J = dets[nJ]

            if I.α == J.α && I.β == J.β
                αindex!(I,αind) 
                βindex!(I,βind) 

                for i=αind
                    γ[i,i,nI,nI] += 1
                end

                for i=βind
                    γ[i,i,nI,nI] += 1
                end
            end

            αexc = αexcitation_level(I,J)
            if αexc > 1
                continue
            end

            βexc = βexcitation_level(I,J)
            if αexc + βexc > 1
                continue
            end

            if αexc == 1
                i, = αexclusive(I,J)
                j, = αexclusive(J,I)
                p = phase(I,J)
                γ[i,j,nI,nJ] = p

            elseif βexc == 1
                i, = βexclusive(I,J)
                j, = βexclusive(J,I)
                p = phase(I,J)
                γ[i,j,nI,nJ] = p
            end
        end
    end

    return γ
end

function build_H_fullγ(dets, h, V)
    nmo = size(h,1)
    γ = get_1p_coupling_coefficients(dets, nmo)
    δ = [i==j ? 1 : 0 for i = 1:nmo, j = 1:nmo]
    @tensoropt H[I,J] := γ[i,j,I,J]*h[i,j] - 0.5*γ[i,l,I,J]*δ[j,k]*V[i,j,k,l] + 0.5*γ[i,j,I,K]*γ[k,l,K,J]*V[i,j,k,l]
    return H
end