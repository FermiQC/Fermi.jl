using Fermi.ConfigurationInteraction.DetOperations

function RFCI(alg::ACI)
    aoints = IntegralHelper{Float64}()
    rhf = Fermi.HartreeFock.RHF(aoints)
    RFCI(aoints, rhf, alg)
end

function RFCI(aoints::IntegralHelper, rhf::Fermi.HartreeFock.RHF, alg::ACI)

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

    ttotal = @elapsed begin
        output("\n    • Computing FCI with the ACI algorithm.\n\n")
        act_range = (Nfrozen+1):(Nactive+Nfrozen)
        σ = Fermi.Options.get("σ")
        γ = Fermi.Options.get("γ")
        pe = Fermi.Options.get("aci_print_screen")
    
        # Determine precision used to represent determinants
        det_size = 
        if Fermi.Options.get("det_size") == 64
            Int64
        elseif Fermi.Options.get("det_size") == 128
            Int128
        else
            throw(Fermi.InvalidFermiOption("Invalid determinant representation $(Fermi.Options.get("det_size"))"))
        end
    
        # Start reference space as HF
        zeroth = repeat('1', Nfrozen)*repeat('1', Int(act_elec/2))
        D0 = Determinant(zeroth, zeroth; precision=det_size)
        P = [D0]
        Pcoef = [1.0]

        E = Eref - Vnuc
        ΔE = 1.0
        ite = 1
    
        output(repeat("=",50))
        Nα = Int(act_elec/2)
        Nβ = Int(act_elec/2)
        M = nothing
        ϵsum = nothing
        ϵest = nothing
        oldP = nothing
        cflag = false
        while true
            output(" → Iteration {}", ite)
            output("Initial model space (P) size: {}", length(P))
    
            output("Generating First Order Interacting Space P ⇒ F ...")
            t = @elapsed F = get_fois(P, Nα, Nβ, act_range)
            output("FOIS size:                {}", length(F))
            output("FOIS contructed in {:5.5f} s.\n", t)
    
            output("Screening FOIS using 2-D Hamiltonian", σ)
            t = @elapsed Fe = ϵI(F, P, Pcoef, E, hp, eri)
            output("Screen complete in {:5.5} s.", t)
            output("Sorting F space...")
            _t = @elapsed begin
                Fperm = zeros(Int, length(Fe))
                sortperm!(Fperm, Fe, by=abs)
                reverse!(Fperm)
                Fe = Fe[Fperm]
                F = F[Fperm]
            end
            output("Sorted in {:5.5f} s.", _t)
            output("Constructing secondary space (Q) by filtering F...")
            t = @elapsed begin
            ϵest = 0.0
            ϵsum = 0.0
            while true
                if length(Fe) == 0 
                    #then no determinants were deemed important - exit ACI
                    break
                end
                if ϵsum ≤ σ
                    ϵest += Fe[end]
                    ϵsum += abs(Fe[end])
                    pop!(Fe)
                    pop!(F)
                else
                    ϵest -= Fe[end]
                    ϵsum -= abs(Fe[end])
                    break
                end
            end
            end
            output("Size of Q: {}", length(Fe))
            output("Secondary space (Q) built in {:5.5f}", t)
    
            ΔE = -E
            M = vcat(P, F)
            output("\nModel space size (P ∪ Q): {}", length(M))
            output("Solving for the model space wave function...", length(M))
            t = @elapsed E, Pcoef, P = update_model_space(M, hp, eri)
            output("Model space updated in {:5.5f}\n", t)
            ΔE += E
            output(" • Model Space Energy           {:15.10f}", E + Vnuc)
            output(" • Energy Change                {:15.10f}", ΔE)
    
            if oldP == Set(P) 
                cflag = true
                break
            end
            ite += 1
            if ite > 30
                break
            end
            oldP = Set(deepcopy(P))
            P, Pcoef = coarse_grain(P, Pcoef, γ, σ)
    
            output("Final coarse grained model space size is {}", length(P))
            output(repeat("=",50)*"\n")
        end
    end #@elapsed
    
    output(repeat("=",50))
    if cflag
        output("🔥🔥🔥🔥🔥 ACI procedure has converged. 🔥🔥🔥🔥🔥")
    else
        output("😲😲😲😲😲 ACI procedure has failed!!!! 😲😲😲😲😲")
    end
    output("Computation finished in {:5.5} seconds.", ttotal)
    output("Model space size: {}", length(M))
    output("E[ACI:{}]     = {:15.10f}", σ, E + Vnuc)
    output("E[ACI:{}]+PT2 = {:15.10f}", σ, E + Vnuc + ϵest)
    output(repeat("=",51)*"\n\n")
    E = (E+Vnuc)
    
    output("\n • Most important determinants:\n")
    output("Coefficient / Determinant / α-Occupancy / β-Occupancy")
    for i in 1:(min(10,length(P)))
        output("{:15.5f}      {}", Pcoef[i], detstring(P[i], Nfrozen+Nactive))
    end
        
    return RFCI(E, E-Eref, Pcoef, P)
end

@fastmath @inbounds function get_fois(dets::Vector{Determinant{T}}, Nα::Int, Nβ::Int, act_range::UnitRange{Int64})::Vector{Determinant{T}} where T <: Integer

    one = typeof(dets[1].α)(1)
    # Ns must be > 1
    Nαvir = length(act_range)-Nα
    Nβvir = length(act_range)-Nβ

    # Estimate FOIS per det
    lf_per_det = (Nα^2*Nαvir^2 + Nα*Nαvir
                       + Nβ^2*Nβvir^2 + Nβ*Nβvir
                       + Nα*Nαvir*Nβ*Nβvir)
    # Estimated total number of determinants (FOIS per det * ndets)
    lf_crit = Int(round(length(dets)*lf_per_det))
    # Preallocate array to hold dummy dets
    fois = [Determinant(0,0) for i=1:lf_crit]

    # Streaming worker-pool: each task keeps its own local scratch buffers
    # (αocc/βocc/αuno/βuno), pulling determinant indices off a shared queue,
    # instead of indexing per-thread arrays by Threads.threadid() (unsafe
    # under task migration / dynamic scheduling). Writes into `fois` are
    # already safe on their own: each _DI owns a disjoint offset range,
    # independent of which task/thread processes it.
    ntasks = Threads.nthreads()
    requests = Channel{Int}(length(dets))
    for _DI in eachindex(dets)
        put!(requests, _DI)
    end
    close(requests)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
            αocc = zeros(Int,Nα)
            βocc = zeros(Int,Nβ)
            αuno = zeros(Int,Nαvir)
            βuno = zeros(Int,Nβvir)
            for _DI in requests
            d = dets[_DI]
            DI = (_DI-1)*lf_per_det + 1
            ct = 0
            αocc!(d, act_range, αocc)
            βocc!(d, act_range, βocc)
            αvir!(d, act_range, αuno)
            βvir!(d, act_range, βuno)
            # Get αα -> αα excitations
            for i in αocc
                for a in αuno
                    newα = (d.α ⊻ (one<<(i-1))) | (one<<(a-1)) 
                    _det = Determinant(newα, d.β)
                    fois[DI+ct] = _det
                    ct += 1
                    for j in αocc
                        if j ≥ i
                            break
                        end
                        for b in αuno
                            if b ≥ a
                                break
                            end
                            newestα = (newα ⊻ (one<<(j-1))) | (one<<(b-1)) 
                            _det = Determinant(newestα, d.β)
                            fois[DI+ct] = _det
                            ct += 1
                        end
                    end
                end
            end
            # Get ββ -> ββ excitations
            for i in βocc
                for a in βuno
                    newβ = (d.β ⊻ (one<<(i-1))) | (one<<(a-1)) 
                    _det = Determinant(d.α, newβ)
                    fois[DI+ct] = _det
                    ct += 1
                    for j in βocc
                        if j ≥ i
                            break
                        end
                        for b in βuno
                            if b ≥ a
                                break
                            end
                            newestβ = (newβ ⊻ (one<<(j-1))) | (one<<(b-1)) 
                            _det = Determinant(d.α, newestβ)
                            fois[DI+ct] = _det
                            ct += 1
                        end
                    end
                end
            end
            # Get αβ -> αβ excitations
            for i in αocc
                for a in αuno
                    for j in βocc
                        for b in βuno
                            newα = (d.α ⊻ (one<<(i-1))) | (one<<(a-1)) 
                            newβ = (d.β ⊻ (one<<(j-1))) | (one<<(b-1)) 
                            _det = Determinant(newα, newβ)
                            fois[DI+ct] = _det
                            ct += 1
                        end
                    end
                end
            end
            end # for _DI in requests
        end #Threads.@spawn
    end
    fois = filter((x)->x != Determinant(0,0), fois)
    fois = Set(fois)
    setdiff!(fois, dets)
    fois = collect(fois)
    return fois 
end

function ϵI(Fdets, P::Vector{Determinant{D}}, Pcoef::Vector{Float64}, Ep::T, h::Array{T,2}, V::Array{T,4}) where {T <: AbstractFloat, D <: Integer}
    Fe = zeros(length(Fdets))
    N = sum(αlist(P[1]))

    # Streaming worker-pool: each task keeps its own local αind/βind scratch
    # buffers, pulling determinant indices off a shared queue. Writes into Fe
    # are inherently safe (each `i` is an independent, disjoint write), so
    # unlike the other fixes in this file this needs no results-merging step
    # -- it only ever needed local scratch buffers instead of indexing a
    # shared per-thread array by Threads.threadid(). (This loop was also
    # missing its Threads.@spawn call entirely -- it ran serially despite
    # looking parallel; restored here.)
    ntasks = Threads.nthreads()
    requests = Channel{Int}(length(Fdets))
    for i in eachindex(Fdets)
        put!(requests, i)
    end
    close(requests)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
        αind = Array{Int64,1}(undef,N)
        βind = Array{Int64,1}(undef,N)
        for i in requests
        D1 = Fdets[i]
        αindex!(D1, αind)
        βindex!(D1, βind)
        Ei = Hd0(αind, βind, h, V)
        Δ = Ei - Ep
        Vint = 0.0
        for j in eachindex(P)
            D2 = P[j]
            αexc = αexcitation_level(D1,D2)
            βexc = βexcitation_level(D1,D2)
            el = αexc + βexc
            if el > 2
                continue
            elseif el == 2
                Vint += Pcoef[j]*Hd2(D1, D2, V, αexc)
            elseif el == 1
                Vint += Pcoef[j]*Hd1(αind, βind, D1, D2, h, V, αexc)
            end
        end

        @fastmath Fe[i] = Δ/2 - √((Δ^2)/4 + Vint^2)
        end
        end #Threads.@spawn
    end

    return Fe
end

function update_model_space(M::Vector{Determinant{D}}, h::Array{T,2}, V::Array{T,4}; complete=true) where {T <: AbstractFloat, D <: Integer}

    if complete
        M = complete_set(M)
    end
    H = get_sparse_hamiltonian_matrix(M, h, V, Fermi.Options.get("cas_cutoff"))

    output("Diagonalizing Hamiltonian...")
    decomp, history = partialschur(H, nev=1, tol=10^-12, which=ArnoldiMethod.LM())
    λ, ϕ = partialeigen(decomp)

    return λ[1], ϕ[:,1], deepcopy(M)
end

function complete_set(dets::Vector{Determinant{T}}) where T <: Integer

    one = typeof(dets[1].α)(1)

    # Streaming worker-pool: each task keeps its own local determinant list
    # and pulls determinants off a shared queue, instead of indexing a
    # per-thread array by Threads.threadid() (unsafe under task migration /
    # dynamic scheduling, and worse here since it's a concurrent push!).
    ntasks = Threads.nthreads()
    requests = Channel{Determinant{T}}(length(dets))
    for d in dets
        put!(requests, d)
    end
    close(requests)

    results = Channel{Vector{Determinant{T}}}(ntasks)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
        local_newdets = Determinant{T}[]
        for d in requests

        asym = d.α ⊻ d.β
        if asym == 0
            continue
        end

        sym = d.α & d.β

        n = count_ones(asym)
        e = Int(n/2)
        idx = Int[]

        str = vcat(repeat([1],e), repeat([0],e))
        perms = multiset_permutations(str, n)

        i = 1
        while (one<<(i-1)) ≤ asym
            if one<<(i-1) & asym ≠ 0
                push!(idx, i)
            end
            i += 1
        end

        for p in perms
            newα = sym
            newβ = sym
            for (x,i) in zip(p,idx)
                if x == 1
                    newα = newα | (one<<(i-1))
                elseif x == 0
                    newβ = newβ | (one<<(i-1))
                end
            end
            push!(local_newdets, Determinant(newα, newβ))
        end
        end
        put!(results, local_newdets)
        end #Threads.@spawn
    end

    newdets = Determinant{T}[]
    for _ in 1:ntasks
        append!(newdets, take!(results))
    end
    return unique(vcat(dets,newdets))
end

function coarse_grain(dets::Vector{Determinant{D}}, C::Vector{T}, γ::Number, σ::Float64) where {T <: AbstractFloat, D <: Integer}

    output("Coarse graining model space for next iteration")
    # Coarse grain
    Cperm = zeros(Int, length(C))
    sortperm!(Cperm, C, by=i->i^2)
    reverse!(Cperm)
    
    Pcoef = C[Cperm]
    P = dets[Cperm]
    
    while true
        if sum(Pcoef[1:end-1].^2) >= 1-γ*σ
            pop!(Pcoef)
            pop!(P)
        else
            break
        end
    end
    return P, Pcoef
end