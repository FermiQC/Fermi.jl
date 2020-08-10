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
    @output "\n\tDF-MP2 energy is {:>14.10f}\n" δ2

    δ3 = 0.0
    v = 1:nvir
    o = 1:nocc

    @output "\tEvaluating 12 diagrams ..."
    #there are 12 of these bad bois
    δ3 = 0.0
    # Diagram 1
    δ3_1 = 0.0
    @output " 1 "
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
            δ3_1 += 2*ruts*_δ3
        end
    end

    # Diagram 2
    δ3_2 = 0.0
    @output "2 "
    for a in o, d in o
        B1 = Bov[:,a,:]
        B2 = Bov[:,d,:]
        @tensor advv[r,s] := B1[Q,r]*B2[Q,s]
        B1 = Boo[:,a,:]
        B2 = Boo[:,d,:]
        @tensor ooad[c,b] := B1[Q,c]*B2[Q,b]
        Δ1 = [(1/(eps[a]+eps[d]-eps[r+nocc]-eps[s+nocc])) for r in v, s in v]
        for c in o, b in o
            cbad = ooad[c,b]
            B1 = Bov[:,c,:]
            B2 = Bov[:,b,:]
            @tensor vvcb[r,s] := B1[Q,r]*B2[Q,s]
            Δ2 = [(1/(eps[c]+eps[b]-eps[r+nocc]-eps[s+nocc])) for r in v, s in v]
            _δ3 = 0.0
            for r in v, s in v
                _δ3 += advv[r,s]*vvcb[r,s]*Δ1[r,s]*Δ2[r,s]
            end
            δ3_2 += 2*cbad*_δ3
        end
    end

    # Diagram 3
    @output "3 "
    δ3_3 = 0.0
    for r in v, c in o
        B1 = Bov[:,:,r]
        B2 = Bov[:,c,:]
        @tensor ocrv[a,t] := B1[Q,a]*B2[Q,t]
        B1 = Bvv[:,r,:]
        B2 = Boo[:,:,c]
        @tensor rovc[t,a] := B1[Q,a]*B2[Q,t]
        Δ1 = [(1/(eps[a]+eps[c]-eps[r+nocc]-eps[t+nocc])) for a in o, t in v]
        for b in o, s in v
            rbsc = rovc[b,s]
            B1 = Bov[:,:,s]
            B2 = Bov[:,b,:]
            @tensor svob[t,a] := B1[Q,a]*B2[Q,t]
            Δ2 = [(1/(eps[a]+eps[b]-eps[s+nocc]-eps[t+nocc])) for a in o, t in v]
            _δ3 = 0.0
            for a in o, t in v
                _δ3 += ocrv[a,t]*svob[t,a]*Δ1[a,t]*Δ2[a,t]
            end
            δ3_3 += -4*rbsc*_δ3
        end
    end

    # Diagram 4
    # r b 
    #   a s
    #      c t
    #           +
    @output "4 "
    δ3_4 = 0.0
    for r in v, b in o
        B1 = Bov[:,b,r]
        B2 = Bov[:,:,:]
        @tensor borv[c,t] := B1[Q]*B2[Q,c,t]
        B1 = Bvv[:,r,:]
        B2 = Boo[:,:,b]
        @tensor rovb[a,s] := B1[Q,s]*B2[Q,a]
        Δ1 = [(1/(eps[b]+eps[c]-eps[r+nocc]-eps[t+nocc])) for c in o, t in v]
        for a in o, s in v
            rasb = rovb[a,s]
            B1 = Bov[:,a,s]
            B2 = Bov[:,:,:]
            @tensor svao[t,c] := B1[Q]*B2[Q,c,t]
            Δ2 = [(1/(eps[a]+eps[c]-eps[s+nocc]-eps[t+nocc])) for c in o, t in v]
            _δ3 = 0.0
            for c in o, t in v
                _δ3 += borv[c,t]*svao[t,c]*Δ1[c,t]*Δ2[c,t]
            end
            δ3_4 += -4*rasb*_δ3
        end
    end

    # Diagram 5
    # c t -> 
    #   b s -> 
    #     a r -> 
    #           +
    @output "5 "
    δ3_5 = 0.0
    for c in o, t in v
        B1 = Bov[:,:,:]
        B2 = Bov[:,c,t]
        @tensor ocvt[a,r] := B1[Q,a,r]*B2[Q]
        B1 = Bov[:,:,:]
        B2 = Bov[:,c,t]
        @tensor otvc[b,s] := B1[Q,b,s]*B2[Q]
        Δ1 = [(1/(eps[a]+eps[c]-eps[r+nocc]-eps[t+nocc])) for a in o, r in v]
        for b in o, s in v
            btsc = otvc[b,s]
            B1 = Bov[:,b,s]
            B2 = Bov[:,:,:]
            @tensor vsob[r,a] := B1[Q]*B2[Q,a,r]
            Δ2 = [(1/(eps[a]+eps[b]-eps[r+nocc]-eps[s+nocc])) for a in o, r in v]
            _δ3 = 0.0
            for a in o, r in v
                _δ3 += ocvt[a,r]*vsob[r,a]*Δ1[a,r]*Δ2[a,r]
            end
            δ3_5 += 8.0*btsc*_δ3
        end
    end

    # Diagram 6
    # c t -> o v
    #   a s -> o v
    #     b r -> o v
    @output "6 "
    δ3_6 = 0.0
    for c in o, t in v
        B1 = Bov[:,c,:]
        B2 = Bov[:,:,t]
        @tensor covt[b,r] := B1[Q,r]*B2[Q,b]
        B1 = Bov[:,:,:]
        B2 = Bov[:,c,t]
        @tensor otvc[a,s] := B1[Q,a,s]*B2[Q]
        Δ1 = [(1/(eps[c]+eps[b]-eps[r+nocc]-eps[t+nocc])) for b in o, r in v]
        for a in o, s in v
            atsc = otvc[a,s]
            B1 = Bov[:,a,:]
            B2 = Bov[:,:,s]
            @tensor vsao[r,b] := B1[Q,r]*B2[Q,b]
            Δ2 = [(1/(eps[a]+eps[b]-eps[r+nocc]-eps[s+nocc])) for b in o, r in v]
            _δ3 = 0.0
            for b in o, r in v
                _δ3 += covt[b,r]*vsao[r,b]*Δ1[b,r]*Δ2[b,r]
            end
            δ3_6 += 2.0*atsc*_δ3
        end
    end

    # Diagram 7
    # a c -> o o
    #   d b -> o o
    #     r s -> v v
    @output "7 "
    δ3_7 = 0.0
    for a in o, c in o
        B1 = Bov[:,a,:]
        B2 = Bov[:,c,:]
        @tensor acrs[r,s] := B1[Q,r]*B2[Q,s]
        B1 = Boo[:,:,a]
        B2 = Boo[:,:,c]
        @tensor dbac[d,b] := B1[Q,d]*B2[Q,b]
        Δ1 = [(1/(eps[a]+eps[c]-eps[r+nocc]-eps[s+nocc])) for r in v, s in v]
        for d in o, b in o
            _dbac = dbac[d,b]
            B1 = Bov[:,d,:]
            B2 = Bov[:,b,:]
            @tensor srdb[s,r] := B1[Q,s]*B2[Q,r]
            Δ2 = [(1/(eps[d]+eps[b]-eps[r+nocc]-eps[s+nocc])) for r in v, s in v]
            _δ3 = 0.0
            for r in v, s in v
                _δ3 += acrs[r,s]*srdb[s,r]*Δ1[r,s]*Δ2[r,s] 
            end
            δ3_7 += -1.0*_dbac*_δ3
        end
    end

    # Diagram 8
    # t r -> v v
    #   u s -> v v
    #     a b -> o o
    @output "8 "
    δ3_8 = 0.0
    for t in v, r in v
        B1 = Bov[:,:,r]
        B2 = Bov[:,:,t]
        @tensor abrt[a,b] := B1[Q,a]*B2[Q,b]
        B1 = Bvv[:,t,:]
        B2 = Bvv[:,r,:]
        @tensor trus[u,s] := B1[Q,u]*B2[Q,s]
        Δ1 = [(1/(eps[a]+eps[b]-eps[t+nocc]-eps[r+nocc])) for a in o, b in o]
        for u in v, s in v
            _trus = trus[u,s]
            B1 = Bov[:,:,u]
            B2 = Bov[:,:,s]
            @tensor usab[a,b] := B1[Q,a]*B2[Q,b]
            Δ2 = [(1/(eps[a]+eps[b]-eps[u+nocc]-eps[s+nocc])) for a in o, b in o]
            _δ3 = 0.0
            for a in o, b in o
                _δ3 += abrt[a,b]*usab[a,b]*Δ1[a,b]*Δ2[a,b]
            end
            δ3_8 += -1.0*_trus*_δ3
        end
    end

    # Diagram 9
    # b r -> o v
    #   s a -> v o
    #     c t -> o v 
    @output "9 "
    δ3_9 = 0.0
    for b in o, r in v
        B1 = Bov[:,b,r]
        B2 = Bov[:,:,:]
        @tensor bcrt[c,t] := B1[Q]*B2[Q,c,t]
        B1 = Boo[:,:,b]
        B2 = Bvv[:,r,:]
        @tensor arbs[a,s] := B1[Q,a]*B2[Q,s]
        Δ1 = [(1/(eps[c]+eps[b]-eps[r+nocc]-eps[t+nocc])) for c in o, t in v]
        for s in v, a in o
            _arbs = arbs[a,s]
            B1 = Bov[:,a,:]
            B2 = Bov[:,:,s]
            @tensor tsac[t,c] := B1[Q,t]*B2[Q,c]
            Δ2 = [(1/(eps[a]+eps[c]-eps[s+nocc]-eps[t+nocc])) for c in o, t in v]
            _δ3 = 0.0
            for c in o, t in v
                _δ3 += bcrt[c,t]*tsac[t,c]*Δ1[c,t]*Δ2[c,t]
            end
            δ3_9 += 2.0*_arbs*_δ3
        end
    end

    # Diagram 10
    # b r -> o v
    #   a s -> o v
    #     c t -> o v
    @output "10 "
    δ3_10 = 0.0
    for b in o, r in v
        B1 = Bov[:,:,r]
        B2 = Bov[:,b,:]
        @tensor cbrt[c,t] := B1[Q,c]*B2[Q,t]
        B1 = Bvv[:,r,:]
        B2 = Boo[:,:,b]
        @tensor rasb[a,s] := B1[Q,s]*B2[Q,a]
        Δ1 = [(1/(eps[c]+eps[b]-eps[r+nocc]-eps[t+nocc])) for c in o, t in v]
        for a in o, s in v
            _rasb = rasb[a,s]
            B1 = Bov[:,a,s]
            B2 = Bov[:,:,:]
            @tensor stac[t,c] := B1[Q]*B2[Q,c,t]
            Δ2 = [(1/(eps[a]+eps[c]-eps[s+nocc]-eps[t+nocc])) for c in o, t in v]
            _δ3 = 0.0
            for c in o, t in v
                _δ3 += cbrt[c,t]*stac[t,c]*Δ1[c,t]*Δ2[c,t]
            end
            δ3_10 += 2.0*_rasb*_δ3
        end
    end

    # Diagram 11
    # a s -> o v 
    #   c t -> o v
    #     r b -> v o
    @output "11 "
    δ3_11 = 0.0
    for a in o, s in v
        B1 = Bov[:,a,:]
        B2 = Bov[:,:,s]
        @tensor abrs[b,r] := B1[Q,r]*B2[Q,b]
        B1 = Bov[:,a,s]
        B2 = Bov[:,:,:]
        @tensor scat[c,t] := B1[Q]*B2[Q,c,t]
        Δ1 = [(1/(eps[a]+eps[b]-eps[r+nocc]-eps[s+nocc])) for  b in o, r in v]
        for c in o, t in v
            _scat = scat[c,t]
            B1 = Bov[:,:,:]
            B2 = Bov[:,c,t]
            @tensor rtbc[r,b] := B1[Q,b,r]*B2[Q]
            Δ2 = [(1/(eps[c]+eps[b]-eps[r+nocc]-eps[t+nocc])) for b in o, r in v]
            _δ3 = 0.0
            for r in v, b in o
                _δ3 += abrs[b,r]*rtbc[r,b]*Δ1[b,r]*Δ2[b,r]
            end
            δ3_11 += -4.0*_scat*_δ3
        end
    end

    # Diagram 12
    # c t -> o v
    #   a s -> o v
    #     r b -> v o
    @output "12 "
    δ3_12 = 0.0
    for c in o, t in v
        B1 = Bov[:,:,:]
        B2 = Bov[:,c,t]
        @tensor bcrt[b,r] := B1[Q,b,r]*B2[Q]
        B1 = Bov[:,:,:]
        B2 = Bov[:,c,t]
        @tensor atsc[a,s] := B1[Q,a,s]*B2[Q]
        Δ1 = [(1/(eps[b]+eps[c]-eps[t+nocc]-eps[r+nocc])) for b in o, r in v]
        for a in o, s in v
            _atsc = atsc[a,s]
            B1 = Bov[:,a,:]
            B2 = Bov[:,:,s]
            @tensor rsab[r,b] := B1[Q,r]*B2[Q,b]
            Δ2 = [(1/(eps[a]+eps[b]-eps[r+nocc]-eps[s+nocc])) for b in o, r in v]
            _δ3 = 0.0
            for r in v, b in o
                _δ3 += bcrt[b,r]*rsab[r,b]*Δ1[b,r]*Δ2[b,r]
            end
            δ3_12 += -4.0*_atsc*_δ3
        end
    end
    end

    δ3 += δ3_1 + δ3_2 + δ3_3 + δ3_4 + δ3_5 + δ3_6 + δ3_7 + δ3_8 + δ3_9 + δ3_10 + δ3_11 + δ3_12
    @output "\n\tDF-MP3 energy is {:>14.10f}\n" δ3

    @output "\tDF-MP3 done in {:>5.2f} s\n" ttotal
    #@output repeat("-",80)*"\n"
    #RMP2{T}(ΔMP2[],Fermi.MemTensor{T}(zeros(T,0,0,0,0)),Fermi.MemTensor{T}(zeros(T,0,0,0,0)))
end

