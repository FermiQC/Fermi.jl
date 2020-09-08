using LinearAlgebra
using LoopVectorization
function RMP2{T}(refwfn::Fermi.HartreeFock.RHF,alg::Fermi.MollerPlesset.DF) where T <: AbstractFloat
    Fermi.MollerPlesset.print_header()
    ints = refwfn.ints
    nocc    = refwfn.ndocc
    nvir    = refwfn.nvir
    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]
    ints.orbs.frozencore = drop_occ
    ints.orbs.frozenvir = drop_vir
    eps = refwfn.eps
    ΔMP2 = 0.0
    @output "TEST\n"
    @output "\tComputing MP2 with DF algorithm\n\n"
    ttotal = @elapsed begin
    @output "\tComputing and transforming integrals ...\n"
    @output "\tComputing integrals ..."
    Fermi.Integrals.aux_ri!(ints)
    @output "\tBasis: {}\n" ints.bname["primary"] 
    @output "\tDF basis: {}\n\n" ints.bname["aux"]
    t = @elapsed ints["B"]
    @output "\t done in {:5.2f} s\n" t

    ints.cache["B"] = convert(Array{T},ints.cache["B"])
    @output "\tTransforming integrals to MO basis ..."
    t = @elapsed Bov = ints["BOV"]
    @output "\t done in {:5.2f} s\n" t
    
    Bis = [zeros(T,size(Bov[:,1,:])) for i=1:Threads.nthreads()]
    Bjs = [zeros(T,size(Bis[1])) for i =1:Threads.nthreads()]
    BABs = [zeros(T,nvir,nvir) for i=1:Threads.nthreads()]
    BBAs = [zeros(T,nvir,nvir) for i=1:Threads.nthreads()]

    ΔMP2s = zeros(T,Threads.nthreads())
    epsd = eps[drop_occ+1:end] #occ - drop_occ
    epsv = eps[nocc+1:end] #virtuals
    _nvir = nvir-drop_vir
    _nocc = nocc-drop_occ
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
    @output "\t Energy computed in {:5.2f} s\n" t

    # generator function for T2 amplitudes, given ijab indices
    function gen(g::Fermi.GeneratedTensor,i,j,a,b)
        naux,nocc,nvir = size(g.data["Bov"])
        sum = 0.0
        for Q in 1:naux
            sum += g.data["Bov"][Q,i,a]*g.data["Bov"][Q,j,b]
        end
        eps = g.data["eps"]
        sum /= eps[i] + eps[j] - eps[a+nocc] - eps[b+nocc]
        sum
    end

    # data required to compute a T2 amplitude
    data = Dict{String,Any}(
                            "Bov" => Bov,
                            "eps" => eps
                           )

    no = convert(UInt16,nocc)
    nv = convert(UInt16,nvir)

    G = Fermi.GeneratedTensor{T}(gen,data,4,[no,no,nv,nv])

    @output "\n\tDF-MP2 energy is {:>14.10f}\n" ΔMP2
    @output "\tDF-MP2 done in {:>5.2f} s\n" ttotal
    @output repeat("-",80)*"\n"

    RMP2{T}(ΔMP2[],G)
end

