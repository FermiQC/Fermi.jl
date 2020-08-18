using LinearAlgebra
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
    @output "\tComputing MP2 with DF algorithm\n\n"
    ttotal = @elapsed begin
    @output "\tComputing and transforming integrals ... "
    t = @elapsed begin
        @output "Clearing cache and setting new basis ... \n"
        Fermi.Integrals.aux_ri!(ints)
        @output "Computing new integrals ...\n"
        ints["B"]
        @output "Converting integrals to {}\n" T
        ints.cache["B"] = convert(Array{T},ints.cache["B"])
        @output "Transforming integrals to MO basis \n"
        Bov = ints["BOV"]
        println(eltype(Bov))
    end
    @output "\tBasis: {}\n" ints.bname["primary"] 
    @output "\tDF basis: {}\n\n" ints.bname["aux"]
    @output " done in {:>5.2f} s\n" t
    Bis = [zeros(T,size(Bov[:,1,:])) for i=1:Threads.nthreads()]
    Bjs = [zeros(T,size(Bis[1])) for i =1:Threads.nthreads()]
    BABs = [zeros(T,nvir,nvir) for i=1:Threads.nthreads()]
    BBAs = [zeros(T,nvir,nvir) for i=1:Threads.nthreads()]

    ΔMP2s = zeros(T,Threads.nthreads())
    Threads.@threads for i in 1:(nocc-drop_occ)
        Bi = Bis[Threads.threadid()]
        Bj = Bjs[Threads.threadid()]
        @views Bi[:,:] = Bov[:,i,:]
        BAB = BABs[Threads.threadid()]
        BBA = BBAs[Threads.threadid()]
        for j in i:(nocc-drop_occ)
            if i != j
                fac = 2
            else
                fac = 1
            end
            @views Bj[:,:] = Bov[:,j,:]
            @tensoropt begin
                BAB[a,b] = Bi[Q,a]*Bj[Q,b]
            end
            #BBA = transpose(BAB)
            transpose!(BBA,BAB)
            for b in 1:nvir-drop_vir
                for a in 1:nvir-drop_vir
                    iajb = BAB[a,b]
                    ibja = BBA[a,b]
                    ΔMP2s[Threads.threadid()] += fac*iajb*(2*iajb - ibja)/(eps[i+drop_occ] + eps[j+drop_occ] - eps[a+nocc] - eps[b+nocc])
                end
            end
        end
    end
    end
    ΔMP2 = sum(ΔMP2s)

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

