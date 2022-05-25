using Combinatorics

function RFCI(alg::RFCIa)
    aoints = IntegralHelper{Float64}()
    rhf = Fermi.HartreeFock.RHF(aoints)

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

    #for iter = 1:10
    ##    
    #end

    σ =  get_σ1(Is, hp, eri, C0) + get_σ3(Is, eri, C0)
    println(sum(σ .* C0) + mol.Vnuc)
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

                        #println("-"^30)
                        #println("$pl $pk $((-1)^(pl+pk))")
                        #println(bitstring(Kβ)[(sz-Nbas):end])
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

# Given a string I, check if the orbital i is occupied
# NOTE: FIRST ORBITAL INDEX = 0
@inline function isocc(I, i)
    return I & (1 << i) != 0
end