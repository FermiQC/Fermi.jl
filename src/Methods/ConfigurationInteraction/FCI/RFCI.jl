using Combinatorics
using KrylovKit

function RFCI(alg::RFCIa)
    aoints = IntegralHelper{Float64}()
    rhf = Fermi.HartreeFock.RHF(aoints)

    ci_header()

    if typeof(aoints.eri_type) === JKFIT
        aoints = IntegralHelper(eri_type=RIFIT())
    elseif Options.get("precision") == "single"
        aoints = IntegralHelper()
    end
    moints = IntegralHelper(orbitals=rhf.orbitals)
    RFCI(moints, aoints, alg)
end

function get_strings(Nelec, Nvir)

    I0 = "0"^Nvir*"1"^Nelec

    perms = multiset_permutations(I0, Nvir+Nelec)

    Is = zeros(Int, length(perms))
    i = 1
    for p in perms
        Is[i] = parse(Int, join(p), base=2)
        i += 1
    end

    return Is
end

function RFCI(moints, aoints, alg)
    mol = moints.molecule
    Nelec = mol.Nα
    Nbas = aoints.orbitals.basisset.nbas
    Nvir = Nbas - Nelec

    # Get I strings
    Is = get_strings(Nelec, Nvir)
    Ns = length(Is)

    # Get integrals
    hp  = Fermi.Integrals.compute!(moints, aoints, "T")
    hp += Fermi.Integrals.compute!(moints, aoints, "V")
    eri = Fermi.Integrals.compute!(moints, "ERI")

    for l = 1:Nbas
        for j = 1:Nbas
            for k = 1:Nbas
                hp[k,l] -= 0.5*eri[k,j,j,l]
            end
        end
    end

    # Trial vector
    C0 = zeros(Ns, Ns)
    C0[1,1] = 1.0

    #f(x) = begin 
    #    @time a = get_σ1(Is, hp, eri, x) 
    #    @time a += get_σ3(Is, eri, x)
    #end
    #σ =  f(C0)
    #println(sum(σ .* C0) + mol.Vnuc)

    #x = eigsolve(f, C0, 1, :LM)

    #x[1][1] + mol.Vnuc

    @time get_σ3(Is, eri, C0)
    @time begin
        tree = build_tree(Is, Nelec, Nvir)
        alt_get_σ3(Is, tree, eri, C0)
    end

end

function build_tree(Is, no, nv)

    N = length(Is)
    tree = zeros(Int, N, no*nv+1)

    off = 1
    for i in 1:N
        Ia = Is[i]
        #println("A $(bitstring(Ia)[62:end])")

        for j in 1:N
            Ib = Is[j]
            #println("B $(bitstring(Ib)[62:end])")

            if count_ones(Ia ⊻ Ib) < 2
                tree[i, off] = j
                off +=1
            end
        end
        off = 1
    end
    return tree
end

function get_σ1(Is, hp, eri, C)

    σ1 = similar(C)
    σ1 .= 0.0

    sz = sizeof(eltype(Is))*8
    Nbas = size(hp, 1)
    F = zeros(length(Is))

    for Iidx = 1:length(Is)
        Iβ = Is[Iidx]
        #println(bitstring(Iβ)[(sz-Nbas):end])
        F .= 0.0

        # Phase associated with annihilating electron l
        pl = -1 
        # loop through l, occ indexes
        for l = 0:(sz - leading_zeros(Iβ) -1)
            if isocc(Iβ,l)
                pl += 1

                # Create new string with l annihilated aₗ|Iβ⟩
                alIβ = Iβ ⊻ (1 << l) # May recycle the bit shift used in isocc

                # Occupied index l found. Now search for an unocupied index k
                pk = 0
                for k = 0:(Nbas - 1)
                    if isocc(alIβ, k)
                        pk += 1
                    else
                        # unoccupied index k found.
                        # Create new string a†ₖaₗ|Iβ⟩   
                        Kβ = alIβ ⊻ (1 << k)

                        Kβidx = findfirst(x->x==Kβ, Is)
                        if isodd(pk+pl)
                            F[Kβidx] -= hp[k+1, l+1]
                        else
                            F[Kβidx] += hp[k+1, l+1]
                        end

                        # Phase associated with annihilating electron j
                        pj = -1 
                        # loop through j, occ indexes
                        for j = 0:(sz - leading_zeros(Kβ) -1)
                            if isocc(Kβ,j)
                                pj += 1

                                # Create new string with j annihilated aⱼ|Kβ⟩
                                ajKβ = Kβ ⊻ (1 << j) # May recycle the bit shift used in isocc

                                # Occupied index j found. Now search for an unocupied index i
                                pi_ = 0
                                for i = 0:(Nbas - 1)
                                    if isocc(ajKβ, i)
                                        pi_ += 1
                                    else
                                        # unoccupied index k found.
                                        # Create new string a†ₖaₗ|Iβ⟩   
                                        Jβ = ajKβ ⊻ (1 << i)

                                        Jβidx = findfirst(x->x==Jβ, Is)
                                        if isodd(pk+pl+pi_+pj)
                                            F[Jβidx] -= 0.5*eri[i+1, j+1, k+1, l+1]
                                        else
                                            F[Jβidx] += 0.5*eri[i+1, j+1, k+1, l+1]
                                        end
                                    end 
                                end # loop over i
                            end
                        end # loop over j
                    end
                end # loop over k
            end
        end # loop over l
        

        σ1[:,Iidx] .= C*F
    end # loop over Is

    return σ1 + transpose(σ1)
end

function get_σ3(Is, eri, C)

    σ3 = similar(C)
    σ3 .= 0.0

    sz = sizeof(eltype(Is))*8

    Nbas = size(eri, 1)
    for αidx in 1:length(Is)
        Iα = Is[αidx]

        ## Loop through l,k to construct Jα string
        # Phase associated with annihilating electron l
        pl = -1 
        # loop through l, occ indexes
        for l = 0:(sz - leading_zeros(Iα) -1)
            if isocc(Iα,l)
                pl += 1

                # Create new string with l annihilated aₗ|Iα⟩
                alIα = Iα ⊻ (1 << l) # May recycle the bit shift used in isocc

                # Occupied index l found. Now search for an unocupied index k
                pk = 0
                for k = 0:(Nbas - 1)
                    if isocc(alIα, k)
                        pk += 1
                    else
                        # unoccupied index k found.
                        # Create new string a†ₖaₗ|Iα⟩   
                        Jα = alIα ⊻ (1 << k)

                        Jαidx = findfirst(x->x==Jα, Is)
                        for βidx in 1:length(Is)
                            Iβ = Is[βidx]

                            ## Loop through i,j to construct Jβ string
                            # Phase associated with annihilating electron j
                            pj = -1 
                            # loop through j, occ indexes
                            for j = 0:(sz - leading_zeros(Iβ) -1)
                                if isocc(Iβ,j)
                                    pj += 1

                                    # Create new string with l annihilated aₗ|Iα⟩
                                    ajIβ = Iβ ⊻ (1 << j) # May recycle the bit shift used in isocc

                                    # Occupied index jl found. Now search for an unocupied index i
                                    pi_ = 0
                                    for i = 0:(Nbas - 1)
                                        if isocc(ajIβ, i)
                                            pi_ += 1
                                        else
                                            # unoccupied index i found.
                                            # Create new string a†ᵢaⱼ|Iβ⟩   
                                            Jβ = ajIβ ⊻ (1 << i)

                                            Jβidx = findfirst(x->x==Jβ, Is)

                                            if isodd(pi_+pj+pl+pk)
                                                σ3[αidx, βidx] -= eri[i+1,j+1,k+1,l+1]*C[Jαidx, Jβidx]
                                            else
                                                σ3[αidx, βidx] += eri[i+1,j+1,k+1,l+1]*C[Jαidx, Jβidx]
                                            end

                                        end
                                    end # loop over i
                                end
                            end # loop over j
                        end # loop over Iβ
                    end
                end # loop over k
            end
        end # loop over l
    end # loop over Iα

    return σ3
end

function alt_get_σ3(Is, tree, eri, C)

    σ3 = similar(C)
    σ3 .= 0.0

    sz = sizeof(eltype(Is))*8
    Nexc = size(tree, 2) + 1

    Nbas = size(eri, 1)
    for αidx in 1:length(Is)
        Iα = Is[αidx]

        # Loop through single excited strings
        # |Iα⟩ = Ekl|Jα⟩
        for exc in 1:Nexc
            Jidx = tree[αidx, exc]
            Jα = Is[Jidx]

            # find k and l
            αxor = Iα ⊻ Jα
            Jαexc = Jα & αxor
            Iαexc = Iα & αxor

            l = sz - leading_zeros(Jαexc)
            k = sz - leading_zeros(Iαexc)

            # Get phase
            p1 = 0
            x0,xf = minmax(l,k)
            for x in (x0+1):(xf-1)
                if Iα & (1 << (x-1)) !== 0 
                    p1 += 1
                end
            end

            for βidx in 1:length(Is)
                Iβ = Is[βidx]

                # Loop through single excited strings
                # |Iβ⟩ = Eij|Jβ⟩
                for exc in 1:Nexc
                    Jidx = tree[βidx, exc]
                    Jβ = Is[Jidx]

                    # find k and l
                    βxor =  Iβ ⊻ Jβ
                    Jβexc = Jβ & βxor
                    Iβexc = Iβ & βxor

                    j = sz - leading_zeros(Jβexc)
                    i = sz - leading_zeros(Iβexc)

                    # Get phase
                    p2 = 0
                    x0,xf = minmax(i,j)
                    for x in (x0+1):(xf-1)
                        if Iβ & (1 << (x-1)) !== 0 
                            p2 += 1
                        end
                    end

                    if isodd(p1+p2)
                        # Need special case for diagonals lol
                        σ3[αidx, βidx] -= eri[i,j,k,l]*C[Jαidx, Jβidx]
                    else
                        σ3[αidx, βidx] += eri[i,j,k,l]*C[Jαidx, Jβidx]
                    end

                end # loop over ij
            end # loop over Iβ
        end # loop over kl
    end # loop over Iα

    return σ3
end

# Given a string I, check if the orbital i is occupied
# NOTE: FIRST ORBITAL INDEX = 0
@inline function isocc(I, i)
    return I & (1 << i) != 0
end