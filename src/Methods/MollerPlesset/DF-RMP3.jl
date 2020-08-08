function RMP3{T}(refwfn::Fermi.HartreeFock.RHF,alg::Fermi.MollerPlesset.DF) where T <: AbstractFloat
    Fermi.MollerPlesset.print_header()
    ints = refwfn.ints
    nocc    = refwfn.ndocc
    nvir    = refwfn.nvir
    eps = refwfn.eps
    δ2 = 0.0
    @output "\tComputing MP3 with DF algorithm\n\n"
    Fermi.Integrals.aux_ri!(ints)
    @output "\tBasis: {}\n" ints.bname["primary"] 
    @output "\tDF basis: {}\n\n" ints.bname["aux"]
    ttotal = @elapsed begin
    @output "\tComputing and transforming integrals ... "
    t = @elapsed begin
        Bov = ints["BOV"]
        Boo = ints["BOO"]
        Bvv = ints["BVV"]
        Bvo = ints["BVO"]
        #Δ   = [(1/(eps[i] + eps[j] - eps[a] - eps[b])) for i=1:nocc,j=1:nocc,a=nocc+1:nocc+nvir,b=nocc+1:nocc+nvir]
    end
    @assert Bov ≈ permutedims(Bvo,(1,3,2))
    @output " done in {:>5.2f} s\n" t
    Bi = zeros(T,size(Bov[:,1,:]))
    Bj = zeros(T,size(Bi))
    BAB = zeros(T,nvir,nvir)
    for i in 1:nocc
        for j in 1:nocc
            Bi[:,:] = Bov[:,i,:]
            Bj[:,:] = Bov[:,j,:]
            @tensoropt begin
                BAB[a,b] = Bi[Q,a]*Bj[Q,b]
            end
            BBA = transpose(BAB)
            for b in 1:nvir
                for a in 1:nvir
                    iajb = BAB[a,b]
                    ibja = BBA[a,b]
                    δ2 += iajb*(2*iajb - ibja)/(eps[i] + eps[j] - eps[a+nocc] - eps[b+nocc])
                end
            end
        end
    end
    end
    @output "\n\tDF-MP2 energy is {:>14.10f}\n" δ2

    δ3 = 0.0
    v = 1:nvir
    o = 1:nocc

    #there are 12 of these bad bois
    #Diagram 1
    for r in v, u in v
        B1 = Bvv[:,:,r]
        B2 = Bvv[:,u,:]
        @tensor ruvv[t,s] := B1[Q,t]*B2[Q,s]
        B1  = Bov[:,:,r]
        B2  = Bov[:,:,u]
        @tensor ooru[a,b] := B1[Q,a]*B2[Q,b]
        Δ1 = [(1/(eps[i] + eps[j] - eps[r+nocc] - eps[u+nocc])) for i in o,j in o]
        for t in v,s in v
            B1 = Bov[:,:,t]
            B2 = Bov[:,:,s]
            @tensor tsoo[a,b] := B1[Q,a]*B2[Q,b]
            ruts = ruvv[t,s]
            Δ2 = [(1/(eps[i] + eps[j] - eps[t+nocc] - eps[s+nocc])) for i in o,j in o]
            #@tensor δ3 += 2*ooru[a,b]*ruts*tsoo[a,b]*Δ[a,b]
            _δ3 = 0.0
            for a in o, b in o
                _δ3 += ooru[a,b]*tsoo[a,b]*Δ1[a,b]*Δ2[a,b]
            end
            δ3 += 2*ruts*_δ3
        end
    end

    #Diagram 2
    #for a in o, d in o
    #    B1 = Boo[:,a,:]
    #    B2 = Boo[:,:,d]
    #    @tensor adoo[r,s]
    @output "\n\tDF-MP3 energy is {:>14.10f}\n" δ3

    @output "\tDF-MP3 done in {:>5.2f} s\n" ttotal
    #@output repeat("-",80)*"\n"
    #RMP2{T}(ΔMP2[],Fermi.MemTensor{T}(zeros(T,0,0,0,0)),Fermi.MemTensor{T}(zeros(T,0,0,0,0)))
end

