using TensorOperations
using LinearAlgebra

function RMP2{T}(refwfn::RHF, ints::IntegralHelper{T}, Alg::RMP2a) where T <: AbstractFloat
    mp_header()

    # Check frozen core and inactive virtual
    nvir = refwfn.nvir
    ndocc = refwfn.ndocc
    core = Options.get("drop_occ")
    inac = Options.get("drop_vir")

    if core ≥ ndocc
        throw(InvalidFermiOption("invalid number of frozen orbitals ($core) for $ndocc doubly occupied orbitals"))
    end
    if inac ≥ nvir
        throw(InvalidFermiOption("invalid number of inactive virtual orbitals ($inac) for $nvir total virtual orbitals"))
    end

    # Build range for MP2 amplitudes
    # Occupied range
    o = (1+core):ndocc
    # Virtual range
    v = (ndocc+1):(ndocc + nvir - inac)

    # Collect RHF information
    orbs = refwfn.orbitals
    C = T.(orbs.C)
    ϵo = T.(orbs.eps[o])
    ϵv = T.(orbs.eps[v])
    o_size = length(ϵo)
    v_size = length(ϵv)

    output("  Starting MP2 computation")
    output(" Number of frozen orbitals:             {:d}" , core)
    output(" Number of inactive orbitals:           {:d}" , inac)
    output(" Number of correlated electron pairs:   {:d}", o_size)
    output(" Number of correlated virtual orbitals: {:d}", v_size)
    output(" ⇒ Total number of MP2 amplitudes:      {:d}\n\n", o_size^2*v_size^2)

    output(repeat("-",80))

    if Options.get("df")
        output("• Performing AO->MO integral transformation...", ending="")

        # TBLIS is great for integrals transformation
        Options.set("tblis", true)
        t = @elapsed begin
        #AOERI = ints["RIERI"]
        #Co = C[:,o]
        #Cv = C[:,v]
        #@tensoropt (P => 100, μ => 50, ν => 50, i => 10, a => 40) begin 
        #    Bov[P,i,a] :=  AOERI[P,μ, ν]*Co[μ, i]*Cv[ν, a]
        #end
        ## Clean data that is not going to be used again
        #AOERI = nothing
        #Co = nothing
        #Cv = nothing
        #delete!(ints, "RIERI")
        #moint = MOIntegralHelper{T, typeof(refwfn.orbitals)}(refwfn.orbitals, ints.auxri, ndocc, ndocc, nvir,nvir, Dict("BOV"=>Bov), false)
        moint = ao_to_mo!(ints, refwfn.orbitals, "BOV")
        end

        # Currently TBLIS is greedy and will not work well inside a threaded loop
        Options.set("tblis", false)

        Bov = moint["BOV"].data
        #Bov = Bov.data
        output("   Done in {:>5.2f} s\n", t)

        output(" Computing DF-MP2 Energy...", ending="")

        aux_size = size(Bov,1)
        # Pre-allocating arrays for threads
        BABs = [zeros(T, v_size, v_size) for i=1:Threads.nthreads()]
        BBAs = [zeros(T, v_size, v_size) for i=1:Threads.nthreads()]

        # Vector containing the energy contribution computed by each thread
        ΔMP2s = zeros(T,Threads.nthreads())
        t = @elapsed begin
            @sync for i in 1:o_size
                Threads.@spawn begin
                id = Threads.threadid()
                @views Bi = Bov[:,i,:]
                BAB = BABs[id]
                BBA = BBAs[id]
                for j in i:o_size
                    if i != j
                        fac = T(2)
                    else
                        fac = one(T)
                    end
                    @views Bj = Bov[:,j,:]

                    @tensor BAB[a,b] = Bi[Q,a]*Bj[Q,b]
                    transpose!(BBA,BAB)

                    eij = ϵo[i] + ϵo[j]
                    dmp2_1 = zero(T)
                    dmp2_2 = zero(T) 
                    @fastmath for b in 1:v_size
                        @inbounds eij_b = eij - ϵv[b]
                        for a in 1:v_size
                            @inbounds begin 
                                iajb = BAB[a,b]
                                d = iajb/(eij_b - ϵv[a])
                                dmp2_1 += d*iajb
                                dmp2_2 += d*BBA[a,b]
                            end
                        end
                    end
                    dmp2 = T(2)*dmp2_1 - dmp2_2
                    ΔMP2s[id] += fac*dmp2
                end
                end
            end
            Emp2 = sum(ΔMP2s)
        end
    else
        output("• Performing AO->MO integral transformation...", ending="")

        # TBLIS is great for integrals transformation
        Options.set("tblis", true)
        t = @elapsed moint = ao_to_mo!(ints, refwfn.orbitals, "OVOV")
        ovov = moint["OVOV"].data
        output("   Done in {:>5.2f} s\n", t)

        output(" Computing MP2 Energy...", ending="")
        t = @elapsed begin
            Emp2 = zero(T)
            for b in 1:nvir
                for a in 1:nvir
                    for j in 1:ndocc
                        for i in 1:ndocc
                            Emp2 += ovov[i,a,j,b]*(2*ovov[i,a,j,b] - ovov[i,b,j,a]) /
                            (ϵo[i]+ϵo[j]-ϵv[a]-ϵv[b])
                        end
                    end
                end
            end
        end
    end
    output("   Done in {:>5.2f} s\n", t)
    output("   @Final RMP2 Correlation Energy {:>20.12f} Eₕ", Emp2)
    output("   @Final RMP2 Total Energy       {:>20.12f} Eₕ", refwfn.energy+Emp2)
    output(repeat("-",80))

    RMP2{T}(refwfn, Emp2, Emp2+refwfn.energy)
end
