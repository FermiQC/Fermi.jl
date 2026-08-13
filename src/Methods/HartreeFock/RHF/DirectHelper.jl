function coulumb_to_fock!(F, D, J, b1::BasisSet, b2::BasisSet)

    # Create a merged basis set 
    # Note that this basis set has repeated atoms. This is done so that
    # we separate auxiliar basis from regular ones. The first N atoms hold the regular basis functions
    # whereas the final N atoms hold auxiliar fitting functions (for a total of 2N atoms).

    atoms = unique(vcat(b1.atoms, b2.atoms))
    basis = vcat(b1.basis, b2.basis)
    b3 = GaussianBasis.BasisSet("$(b1.name*b2.name)", atoms, basis)
    Nvals = GaussianBasis.num_basis.(b3.basis)

    #b3 = GaussianBasis.BasisSet(b1, b2)
    #Nvals = [GaussianBasis.Libcint.CINTcgtos_spheric(i-1, b3.lc_bas) for i = 1:b3.nshells]
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:b3.nshells]
    Nmax = maximum(Nvals)
    ntasks = Threads.nthreads()

    # Build intermediate c
    # Each worker owns a private `buf` (not indexed by Threads.threadid(), which is
    # unsafe if a task migrates threads mid-execution) and pulls shell indices off a
    # shared channel. Different A values write to disjoint c[Aoff+1:Aoff+Nα] ranges,
    # so no synchronization is needed for the accumulation itself.
    c = zeros(b2.nbas)
    Arange = Channel{Int}(b3.nshells - b1.nshells)
    for A in (b1.nshells+1):b3.nshells # A runs through Auxiliar shells only
        put!(Arange, A)
    end
    close(Arange)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
            buf = zeros(Cdouble, Nmax^3)
            for A in Arange
                Nα = Nvals[A]                  # Number of basis in the A shell
                Aoff = ao_offset[A] - b1.nbas  # Number of basis set before the A shell
                for μ in 1:b1.nshells          # Shell index μ
                    Nμ = Nvals[μ]              # Number of basis in the μ shell
                    μoff = ao_offset[μ]        # Number of basis set before the μ shell
                    for ν in μ:b1.nshells      # Shell index ν
                        Nν = Nvals[ν]          # Number of basis in the ν shell
                        νoff = ao_offset[ν]    # Number of basis set before the ν shell
                        f = μ != ν ? 2.0 : 1.0

                        # Compute integral for the (μ,ν,A) triplet shell
                        GaussianBasis.ERI_2e3c!(buf, b3, μ, ν, A)

                        for Am in 1:Nα
                            for νm in 1:Nν
                                _ν = νoff + νm  # Basis set index ν
                                for μm in 1:Nμ
                                    _μ = μoff + μm      # Basis set index μ
                                    γμν = f*D[_μ,_ν]
                                    # C(A) = γ(μν)*(μν|A)
                                    c[Aoff+Am] += γμν*buf[μm + Nμ*(νm-1) + Nμ*Nν*(Am-1)]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Build intermediate d
    # d(A) = J(AB)*C(B)
    d = J*c

    # Add Coulumb contribution to F
    # Same worker-pool treatment. Different μ values write to disjoint (row,col) pairs
    # of F (each F[i,j] is owned by the task whose μ-shell equals min(shell(i),shell(j))),
    # so this was already race-free across tasks -- only the shared `buf` was unsafe.
    μrange = Channel{Int}(b1.nshells)
    for μ in 1:b1.nshells
        put!(μrange, μ)
    end
    close(μrange)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
            buf = zeros(Cdouble, Nmax^3)
            for μ in μrange
                Nμ = Nvals[μ]              # Number of basis in the μ shell
                μoff = ao_offset[μ]        # Number of basis set before the μ shell
                for ν in μ:b1.nshells      # Shell index ν
                        Nν = Nvals[ν]          # Number of basis in the ν shell
                        νoff = ao_offset[ν]    # Number of basis set before the ν shell
                    for A in (b1.nshells+1):b3.nshells # A runs through Auxiliar shells only
                        Nα = Nvals[A]                  # Number of basis in the A shell
                        Aoff = ao_offset[A] - b1.nbas  # Number of basis set before the A shell
                        # Compute integral for the (μ,ν,A) triplet shell
                        GaussianBasis.ERI_2e3c!(buf, b3, μ, ν, A)

                        for Am in 1:Nα
                            for νm in 1:Nν
                                _ν = νoff + νm  # Basis set index ν
                                for μm in 1:Nμ
                                    _μ = μoff + μm      # Basis set index μ
                                    if _μ > _ν
                                        continue
                                    end
                                    # F += 2J  where J = d(A)*(μν|A)
                                    Fμν = 2.0*d[Aoff+Am]*buf[μm + Nμ*(νm-1) + Nμ*Nν*(Am-1)]
                                    F[_μ, _ν] += Fμν
                                    if _μ != _ν
                                        F[_ν, _μ] += Fμν
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function exchange_to_fock!(F, C, J, ndocc, b1::BasisSet, b2::BasisSet)

    # Create a merged basis set 
    # Note that this basis set has repeated atoms. This is done so that
    # we separate auxiliar basis from regular ones. The first N atoms hold the regular basis functions
    # whereas the final N atoms hold auxiliar fitting functions (for a total of 2N atoms).

    atoms = unique(vcat(b1.atoms, b2.atoms))
    basis = vcat(b1.basis, b2.basis)
    b3 = GaussianBasis.BasisSet("$(b1.name*b2.name)", atoms, basis)
    Nvals = GaussianBasis.num_basis.(b3.basis)
    ao_offset = [sum(Nvals[1:(i-1)]) for i = 1:b3.nshells]
    Nmax = maximum(Nvals)
    ntasks = Threads.nthreads()

    # Build half-transformed (μi|A)
    # Same worker-pool treatment as coulumb_to_fock! -- private per-task `buf`,
    # different A values write to disjoint W[:, :, _A] slices so no synchronization
    # is needed for the accumulation itself.
    W = zeros(b1.nbas, ndocc, b2.nbas)
    Arange = Channel{Int}(b3.nshells - b1.nshells)
    for A in (b1.nshells+1):b3.nshells # A runs through Auxiliar shells only
        put!(Arange, A)
    end
    close(Arange)

    @sync for _ in 1:ntasks
        Threads.@spawn begin
            buf = zeros(Cdouble, Nmax^3)
            for A in Arange
                Nα = Nvals[A]                  # Number of basis in the A shell
                Aoff = ao_offset[A] - b1.nbas  # Number of basis set before the A shell
                for μ in 1:b1.nshells          # Shell index μ
                    Nμ = Nvals[μ]              # Number of basis in the μ shell
                    μoff = ao_offset[μ]        # Number of basis set before the μ shell
                    for ν in μ:b1.nshells      # Shell index ν
                        Nν = Nvals[ν]          # Number of basis in the ν shell
                        νoff = ao_offset[ν]    # Number of basis set before the ν shell

                        # Compute integral for the (μ,ν,A) triplet shell
                        GaussianBasis.ERI_2e3c!(buf, b3, μ, ν, A)

                        for Am in 1:Nα
                            _A = Aoff + Am
                            for νm in 1:Nν
                                _ν = νoff + νm      # Basis set index ν
                                for μm in 1:Nμ
                                    _μ = μoff + μm  # Basis set index μ
                                    if _μ > _ν
                                        break
                                    end
                                    for i = 1:ndocc
                                        W[_μ, i, _A] += C[_ν,i]*buf[μm + Nμ*(νm-1) + Nμ*Nν*(Am-1)]
                                        if _μ != _ν
                                            W[_ν, i, _A] += C[_μ,i]*buf[μm + Nμ*(νm-1) + Nμ*Nν*(Am-1)]
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end #spawn
    end #sync

    @tensoropt F[μ,ν] -=  W[μ,i,A]*J[A,B]*W[ν,i,B]
end