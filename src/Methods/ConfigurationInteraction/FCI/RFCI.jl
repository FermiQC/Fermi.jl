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
end

function σ1(Is, hp, eri)

    sz = sizeof(eltype(Is))*8
    Nbas = size(hp, 1)
    F = zeros(length(Is))
    for Iβ = Is
        #println(bitstring(Iβ)[(sz-Nbas):end])
        F .= 0.0

        lz = leading_zeros(Iβ)
        # Phase associated with annihilating electron l
        pl = -1 
        # loop through l, occ indexes
        for l = 0:(sz - lz -1)
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

                        #println("-"^30)
                        #println("$pl $pk $((-1)^(pl+pk))")
                        #println(bitstring(Kβ)[(sz-Nbas):end])
                    end
                end
            end
        end
    end

        


end

# Given a string I, check if the orbital i is occupied
# NOTE: FIRST ORBITAL INDEX = 0
@inline function isocc(I, i)
    return I & (1 << i) != 0
end