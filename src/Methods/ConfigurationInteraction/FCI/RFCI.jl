using LinearAlgebra
using Combinatorics
using KrylovKit
using Molecules

# Given a string I, check if the orbital i is occupied
# NOTE: FIRST ORBITAL INDEX = 0
@inline function isocc(I, i)
    return I & (1 << i) != 0
end

function detstring(I, n = 7)
    return reverse(bitstring(I))[1:n]
end

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

function get_strings(Nelec, Nvir, frozen)

    core = "1"^frozen
    I0 = "0"^Nvir*"1"^Nelec

    perms = multiset_permutations(I0, Nvir+Nelec)

    Is = zeros(Int, length(perms))
    i = 1
    for p in perms
        Is[i] = parse(Int, join(p)*core, base=2)
        i += 1
    end

    return Is
end

function RFCI(moints, aoints, alg::RFCIa)
    Eref = moints.orbitals.sd_energy
    mol = moints.molecule
    Nelec = mol.Nα
    Nbas = aoints.orbitals.basisset.nbas
    Nvir = Nbas - Nelec

    Nfrozen = Options.get("drop_occ")
    Ninac = Options.get("drop_vir")
    Nactive = Nbas - Nfrozen - Ninac

    # Get I strings
    Is = get_strings(Nelec - Nfrozen, Nvir - Ninac, Nfrozen)
    Ns = length(Is)

    output(" => Active Space Information ({:d}e, {:d}o)", Nelec-Nfrozen, Nactive)
    output(" • Number of Electrons:      {:>5d}", Nelec)
    output(" • Active Electrons:         {:>5d}", Nelec - Nfrozen)
    output(" • Number of Orbitals:       {:>5d}", Nbas)
    output(" • Active Orbitals:          {:>5d}", Nactive)
    output("\n => CI dimension")
    output(" • Number of Strings:      {:>5d}", Ns)
    output(" • Number of Determinants: {:>5d}", Ns^2)

    # Get integrals
    output("Transforming integrals... ", ending="")
    t = @elapsed begin
        hp  = Fermi.Integrals.compute!(moints, aoints, "T")
        hp += Fermi.Integrals.compute!(moints, aoints, "V")
        eri = Fermi.Integrals.compute!(moints, "ERI")
    end

    for l = 1:Nbas
        for j = 1:(Nbas-Ninac)
            for k = 1:Nbas
                hp[k,l] -= 0.5*eri[k,j,j,l]
            end
        end
    end
    output("Done in {:10.5f} seconds.", t)

    # Trial vector
    C0 = zeros(Ns, Ns)
    C0[1,1] = 1.0

    tree = build_tree(Is, Nelec - Nfrozen, Nvir - Ninac)
    
    linmap(x) = begin 
        a = get_σ1(Is, hp, eri, x, tree, Nelec)
        a += get_σ3(Is, tree, eri, x, Nelec)
    end
    output("\nStarting eigsolve routine\n")
    klv = eigsolve(linmap, C0, 1, :LM; issymmetric=true, tol=1e-8)

    output("\nKrylov Solver Summary:\n {}", string(klv[3]))

    Efci = klv[1][1] + Molecules.nuclear_repulsion(mol.atoms)
    output("Final FCI Energy: {:15.10f}", Efci)

    return RFCI(Efci, Efci - Eref)
end

function build_tree(Is, no, nv)

    N = length(Is)
    tree = zeros(Int, N, no*nv+1)

    off = 1
    for i in 1:N
        Ia = Is[i]

        for j in 1:N
            Ib = Is[j]

            if count_ones(Ia ⊻ Ib) ≤ 2
                tree[i, off] = j
                off +=1
            end
        end
        off = 1
    end
    return tree
end

function get_σ1(Is, hp, eri, C, tree, Nelec)

    σ1 = similar(C)
    σ1 .= 0.0

    Nexc = size(tree, 2)
    sz = sizeof(eltype(Is))*8
    F = zeros(length(Is))
    Iocc_ind = zeros(Int, Nelec)
    Kocc_ind = zeros(Int, Nelec)

    # Loop over Iβ
    for Iidx in 1:length(Is)

        F .= 0.0      # Set array F(Jβ) = 0
        I = Is[Iidx]

        # Loop over excitations Ekl from Iβ
        for exc in 1:Nexc

            # Get K = Ekl|Iβ⟩
            Kidx = tree[Iidx, exc]
            K = Is[Kidx]

            # find k and l
            IKxor = I ⊻ K

            # If I == K
            if IKxor == 0
                n = 1
                for k = 1:(sz - leading_zeros(I))
                    if isocc(I, k-1)
                        F[Kidx] += hp[k,k]
                        Iocc_ind[n] = k
                        n += 1
                    end # if it's occ
                end # look k vals
            else
                Kexc = K & IKxor
                Iexc = I & IKxor

                l = sz - leading_zeros(Iexc)
                k = sz - leading_zeros(Kexc)
                
                # Get phase
                p1 = 0
                x0,xf = minmax(l,k)
                for x in (x0+1):(xf-1)
                    if I & (1 << (x-1)) !== 0 
                        p1 += 1
                    end
                end

                if isodd(p1)
                    F[Kidx] -= hp[k,l]
                else
                    F[Kidx] += hp[k,l]
                end
            end # I = K

            # Loop over excitations Eij from K
            for exc2 in 1:Nexc

                # Get J = Ekl|K⟩
                Jidx = tree[Kidx, exc2]
                J = Is[Jidx]

                # find i and j
                KJxor = K ⊻ J

                # If J == K
                if KJxor == 0
                    n = 1
                    for k = 1:(sz - leading_zeros(K))
                        if isocc(K, k-1)
                            Kocc_ind[n] = k
                            n += 1
                        end # if it's occ
                    end # look k vals
                else
                    Jexc = J & KJxor
                    Kexc = K & KJxor

                    j = sz - leading_zeros(Kexc)
                    i = sz - leading_zeros(Jexc)

                    # Get phase
                    p2 = 0
                    x0,xf = minmax(i,j)
                    for x in (x0+1):(xf-1)
                        if K & (1 << (x-1)) !== 0 
                            p2 += 1
                        end
                    end
                end # J == K

                # F(J) += 0.5 * sign *[ij|kl]
                if I != K != J
                    if isodd(p1+p2)
                        F[Jidx] -= 0.5*eri[i,j,k,l]
                    else
                        F[Jidx] += 0.5*eri[i,j,k,l]
                    end
                elseif K != J
                    if isodd(p2)
                        for k = Iocc_ind
                            F[Jidx] -= 0.5*eri[i,j,k,k]
                        end
                    else
                        for k = Iocc_ind
                            F[Jidx] += 0.5*eri[i,j,k,k]
                        end
                    end
                elseif I != K
                    if isodd(p1)
                        for i = Kocc_ind
                            F[Jidx] -= 0.5*eri[i,i,k,l]
                        end
                    else
                        for i = Kocc_ind
                            F[Jidx] += 0.5*eri[i,i,k,l]
                        end
                    end
                else
                    for i = Iocc_ind
                        for k = Kocc_ind
                            F[Jidx] += 0.5*eri[i,i,k,k]
                        end
                    end
                end # F(J) += 0.5 * sign *[ij|kl]
            end # Loop over excitations Eij from K
        end # Loop over excitations Ekl from Iβ
        σ1[:,Iidx] .= C*F
    end # Loop over I 

    return σ1 + transpose(σ1)
end


function get_σ3(Is, tree, eri, C, Nelec)

    σ3 = similar(C)
    σ3 .= 0.0

    sz = sizeof(eltype(Is))*8
    Nexc = size(tree, 2)

    αocc_ind = zeros(Int, Nelec)
    βocc_ind = zeros(Int, Nelec)

    for αidx in 1:length(Is)
        Iα = Is[αidx]

        # Loop through single excited strings
        # |Iα⟩ = Ekl|Jα⟩
        for exc in 1:Nexc
            Jαidx = tree[αidx, exc]
            Jα = Is[Jαidx]

            # find k and l
            αxor = Iα ⊻ Jα

            # If Iα = Jα
            if αxor == 0
                n = 1
                for k = 0:(sz - leading_zeros(Iα) -1)
                    if isocc(Iα, k)
                        αocc_ind[n] = k+1
                        n += 1
                    end # if it's occ
                end # look k vals
            end # Iα = Jα


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
                    Jβidx = tree[βidx, exc]
                    Jβ = Is[Jβidx]

                    # find k and l
                    βxor =  Iβ ⊻ Jβ

                    # If Iβ = Jβ
                    if βxor == 0
                        n = 1
                        for i = 0:(sz - leading_zeros(Iβ) -1)
                            if isocc(Iβ, i)
                                βocc_ind[n] = i+1
                                n += 1 
                            end # if it's occ
                        end # look i vals
                    end # Iβ = Jβ

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

                    if (i != j) & (k != l)
                        if isodd(p1+p2)
                            σ3[αidx, βidx] -= eri[i,j,k,l]*C[Jαidx, Jβidx]
                        else
                            σ3[αidx, βidx] += eri[i,j,k,l]*C[Jαidx, Jβidx]
                        end
                    elseif (i != j)
                        if isodd(p2)
                            for k = αocc_ind
                                σ3[αidx, βidx] -= eri[i,j,k,k]*C[Jαidx, Jβidx]
                            end
                        else
                            for k = αocc_ind
                                σ3[αidx, βidx] += eri[i,j,k,k]*C[Jαidx, Jβidx]
                            end
                        end
                    elseif (k != l)
                        if isodd(p1)
                            for i = βocc_ind
                                σ3[αidx, βidx] -= eri[i,i,k,l]*C[Jαidx, Jβidx]
                            end
                        else
                            for i = βocc_ind
                                σ3[αidx, βidx] += eri[i,i,k,l]*C[Jαidx, Jβidx]
                            end
                        end
                    else
                        for i = βocc_ind
                            for k = αocc_ind
                                σ3[αidx, βidx] += eri[i,i,k,k]*C[Jαidx, Jβidx]
                            end
                        end
                    end
                end # loop over ij
            end # loop over Iβ
        end # loop over kl
    end # loop over Iα

    return σ3
end

#function alt_get_σ1(Is, hp, eri, C, Nfrozen, Ninac)
#
#    σ1 = similar(C)
#    σ1 .= 0.0
#
#    sz = sizeof(eltype(Is))*8
#    Nbas = size(hp, 1) - Ninac
#    F = zeros(length(Is))
#
#    for Iidx = 1:length(Is)
#        Iβ = Is[Iidx]
#        F .= 0.0
#
#        # Phase associated with annihilating electron l
#        pl = -1 
#        # loop through l, occ indexes
#        for l = (0+Nfrozen):(sz - leading_zeros(Iβ) -1)
#            if isocc(Iβ,l)
#                pl += 1
#
#                # Create new string with l annihilated aₗ|Iβ⟩
#                alIβ = Iβ ⊻ (1 << l) # May recycle the bit shift used in isocc
#
#                # Occupied index l found. Now search for an unocupied index k
#                pk = 0
#                for k = (0+Nfrozen):(Nbas - 1)
#                    if isocc(alIβ, k)
#                        pk += 1
#                    else
#                        # unoccupied index k found.
#                        # Create new string a†ₖaₗ|Iβ⟩   
#                        Kβ = alIβ ⊻ (1 << k)
#
#                        Kβidx = findfirst(x->x==Kβ, Is)
#                        if isodd(pk+pl)
#                            F[Kβidx] -= hp[k+1, l+1]
#                        else
#                            F[Kβidx] += hp[k+1, l+1]
#                        end
#
#                        # Phase associated with annihilating electron j
#                        pj = -1 
#                        # loop through j, occ indexes
#                        for j = (0+Nfrozen):(sz - leading_zeros(Kβ) -1)
#                            if isocc(Kβ,j)
#                                pj += 1
#
#                                # Create new string with j annihilated aⱼ|Kβ⟩
#                                ajKβ = Kβ ⊻ (1 << j) # May recycle the bit shift used in isocc
#
#                                # Occupied index j found. Now search for an unocupied index i
#                                pi_ = 0
#                                for i = (0+Nfrozen):(Nbas - 1)
#                                    if isocc(ajKβ, i)
#                                        pi_ += 1
#                                    else
#                                        # unoccupied index k found.
#                                        # Create new string a†ₖaₗ|Iβ⟩   
#                                        Jβ = ajKβ ⊻ (1 << i)
#
#                                        Jβidx = findfirst(x->x==Jβ, Is)
#                                        if isodd(pk+pl+pi_+pj)
#                                            F[Jβidx] -= 0.5*eri[i+1, j+1, k+1, l+1]
#                                        else
#                                            F[Jβidx] += 0.5*eri[i+1, j+1, k+1, l+1]
#                                        end
#                                    end 
#                                end # loop over i
#                            end
#                        end # loop over j
#                    end
#                end # loop over k
#            end
#        end # loop over l
#        
#
#        σ1[:,Iidx] .= C*F
#    end # loop over Is
#
#    return σ1 + transpose(σ1)
#end
#
#function alt_get_σ3(Is, eri, C)
#
#    σ3 = similar(C)
#    σ3 .= 0.0
#
#    sz = sizeof(eltype(Is))*8
#
#    Nbas = size(eri, 1)
#    for αidx in 1:length(Is)
#        Iα = Is[αidx]
#
#        ## Loop through l,k to construct Jα string
#        # Phase associated with annihilating electron l
#        pl = -1 
#        # loop through l, occ indexes
#        for l = 0:(sz - leading_zeros(Iα) -1)
#            if isocc(Iα,l)
#                pl += 1
#
#                # Create new string with l annihilated aₗ|Iα⟩
#                alIα = Iα ⊻ (1 << l) # May recycle the bit shift used in isocc
#
#                # Occupied index l found. Now search for an unocupied index k
#                pk = 0
#                for k = 0:(Nbas - 1)
#                    if isocc(alIα, k)
#                        pk += 1
#                    else
#                        # unoccupied index k found.
#                        # Create new string a†ₖaₗ|Iα⟩   
#                        Jα = alIα ⊻ (1 << k)
#
#                        Jαidx = findfirst(x->x==Jα, Is)
#                        for βidx in 1:length(Is)
#                            Iβ = Is[βidx]
#
#                            ## Loop through i,j to construct Jβ string
#                            # Phase associated with annihilating electron j
#                            pj = -1 
#                            # loop through j, occ indexes
#                            for j = 0:(sz - leading_zeros(Iβ) -1)
#                                if isocc(Iβ,j)
#                                    pj += 1
#
#                                    # Create new string with l annihilated aₗ|Iα⟩
#                                    ajIβ = Iβ ⊻ (1 << j) # May recycle the bit shift used in isocc
#
#                                    # Occupied index jl found. Now search for an unocupied index i
#                                    pi_ = 0
#                                    for i = 0:(Nbas - 1)
#                                        if isocc(ajIβ, i)
#                                            pi_ += 1
#                                        else
#                                            # unoccupied index i found.
#                                            # Create new string a†ᵢaⱼ|Iβ⟩   
#                                            Jβ = ajIβ ⊻ (1 << i)
#
#                                            Jβidx = findfirst(x->x==Jβ, Is)
#
#                                            if isodd(pi_+pj+pl+pk)
#                                                σ3[αidx, βidx] -= eri[i+1,j+1,k+1,l+1]*C[Jαidx, Jβidx]
#                                            else
#                                                σ3[αidx, βidx] += eri[i+1,j+1,k+1,l+1]*C[Jαidx, Jβidx]
#                                            end
#
#                                        end
#                                    end # loop over i
#                                end
#                            end # loop over j
#                        end # loop over Iβ
#                    end
#                end # loop over k
#            end
#        end # loop over l
#    end # loop over Iα
#
#    return σ3
#end

