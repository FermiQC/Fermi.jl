using LinearAlgebra
using LoopVectorization

function RMP2{T}(refwfn::RHF, ints::IntegralHelper{T}) where T <: AbstractFloat
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

    output("  Starting MP2 computation")
    output(" Number of frozen orbitals:             {:d}" , core)
    output(" Number of inactive orbitals:           {:d}" , inac)
    output(" Number of correlated electron pairs:   {:d}", length(o))
    output(" Number of correlated virtual orbitals: {:d}", length(v))
    output(" ⇒ Total number of MP2 amplitudes:      {:d}\n\n", length(o)^2*length(v)^2)

    output(repeat("-",80))
    output("\tTransforming integrals to MO basis ...")
    t = @elapsed Bov = Fermi.Integrals.ao_to_mo_rieri!(ints, C[:,o], C[:,v])
    output("\t done in {:5.2f} s\n", t)
    
    Bis = [zeros(T,size(Bov[:,1,:])) for i=1:Threads.nthreads()]
    Bjs = [zeros(T,size(Bis[1])) for i =1:Threads.nthreads()]
    BABs = [zeros(T,nvir,nvir) for i=1:Threads.nthreads()]
    BBAs = [zeros(T,nvir,nvir) for i=1:Threads.nthreads()]

    ΔMP2s = zeros(T,Threads.nthreads())
    epsd = eps[core+1:end] #occ - drop_occ
    epsv = eps[ndocc+1:end] #virtuals
    _nvir = nvir-inac
    _nocc = ndocc-core
    t = @elapsed begin
        @sync for i in 1:_nocc
            Threads.@spawn begin
            id = Threads.threadid()
            Bi = Bis[id]
            Bj = Bjs[id]
            @views Bi[:,:] = Bov[:,i,:]
            BAB = BABs[id]
            BBA = BBAs[id]
            for j in i:_nocc
                if i != j
                    fac = 2
                else
                    fac = 1
                end
                @views Bj[:,:] .= Bov[:,j,:]
                @tensoropt begin
                    BAB[a,b] = Bi[Q,a]*Bj[Q,b]
                end
                transpose!(BBA,BAB)
                eij = epsd[i] + epsd[j]
                dmp2_1 = 0.0
                dmp2_2 = 0.0
                @fastmath for b in 1:_nvir
                    @inbounds eij_b = eij - epsv[b]
                    for a in 1:_nvir
                        @inbounds iajb = BAB[a,b]
                        #@inbounds ibja = BBA[a,b]
                        d = iajb/(eij_b - epsv[a])
                        @inbounds dmp2_1 += d*iajb
                        @inbounds dmp2_2 += d*BBA[a,b]
                    end
                end
                dmp2 = 2*dmp2_1 - dmp2_2
                ΔMP2s[id] += fac*dmp2
            end
            end
        end
        end
        ΔMP2 = sum(ΔMP2s)
    end
    output("\t Energy computed in {:5.2f} s\n", t)

    output("\n\tDF-MP2 energy is {:>14.10f}\n", ΔMP2)
    output("\tDF-MP2 done in {:>5.2f} s\n", ttotal)
    output(repeat("-",80)*"\n")
end

