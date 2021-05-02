function RCCSDpT(Alg::RpTa)
    aoints = IntegralHelper()
    rhf = RHF(aoints)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    ccsd = RCCSD(moints, aoints)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(ccsd::RCCSD, moints::IntegralHelper{T,E,O}, Alg::RpTa) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    output("\n   • Perturbative Triples Started\n")

    T1 = ccsd.T1
    T2 = ccsd.T2

    Vovvv = moints["OVVV"]
    Vooov = moints["OOOV"]
    Vovov = moints["OVOV"]

    o,v = size(T1)

    fo = moints["Fii"]
    fv = moints["Faa"]

    # Pre-allocate Intermediate arrays
    Ws  = [Array{T}(undef, v,v,v) for i = 1:Threads.nthreads()] 
    Vs  = [Array{T}(undef, v,v,v) for i = 1:Threads.nthreads()] 
    Evals = zeros(T, Threads.nthreads())

    output("Computing energy contribution from occupied orbitals:")
    @sync for i in 1:o
        Threads.@spawn begin
            output("→ Orbital {} of {}", i, o)
            Vovvv_1i = view(Vovvv, i,:,:,:)

            T2_1i    = view(T2, i,:,:,:)
            T1_1i    = view(T1, i, :)
            for j in 1:i
                Dij = fo[i] + fo[j]

                Vovvv_1j    = view(Vovvv, j,:,:,:)
                Vovov_1i_3j = view(Vovov, i,:,j,:)

                Vooov_2i_3j = view(Vooov,:,i,j,:)
                Vooov_2j_3i = view(Vooov,:,j,i,:)

                T2_1i_2j    = view(T2, i,j,:,:)
                T2_1j_2i    = view(T2, j,i,:,:)
                T2_1j       = view(T2, j,:,:,:)
                T1_1j       = view(T1, j, :)
                δij = i == j
                for k in 1:j
                    Dijk = Dij + fo[k]

                    Vovvv_1k    = view(Vovvv, k,:,:,:)

                    Vooov_2j_3k = view(Vooov,:,j,k,:)
                    Vooov_2k_3j = view(Vooov,:,k,j,:)
                    Vooov_2k_3i = view(Vooov,:,k,i,:)
                    Vooov_2i_3k = view(Vooov,:,i,k,:)

                    T2_1k_2j = view(T2, k,j,:,:)
                    T2_1k_2i = view(T2, k,i,:,:)
                    T2_1i_2k = view(T2, i,k,:,:)
                    T2_1j_2k = view(T2, j,k,:,:)

                    T2_1k = view(T2, k,:,:,:)
                    T1_1k = view(T1, k, :)

                    Vovov_1j_3k = view(Vovov, j,:,k,:)
                    Vovov_1i_3k = view(Vovov, i,:,k,:)

                    # Build intermediates W and V for the set i j k
                    W = Ws[Threads.threadid()]
                    V = Vs[Threads.threadid()]
                    @tensoropt begin
                        W[a,b,c]  =  Vovvv_1i[a,b,d]*T2_1k_2j[c,d] - Vooov_2j_3k[l,c]*T2_1i[l,a,b]  # ijk abc
                        W[a,b,c] +=  Vovvv_1i[a,c,d]*T2_1j_2k[b,d] - Vooov_2k_3j[l,b]*T2_1i[l,a,c]  # ikj acb
                        W[a,b,c] +=  Vovvv_1k[c,a,d]*T2_1j_2i[b,d] - Vooov_2i_3j[l,b]*T2_1k[l,c,a]  # kij cab
                        W[a,b,c] +=  Vovvv_1k[c,b,d]*T2_1i_2j[a,d] - Vooov_2j_3i[l,a]*T2_1k[l,c,b]  # kji cba
                        W[a,b,c] +=  Vovvv_1j[b,c,d]*T2_1i_2k[a,d] - Vooov_2k_3i[l,a]*T2_1j[l,b,c]  # jki bca
                        W[a,b,c] +=  Vovvv_1j[b,a,d]*T2_1k_2i[c,d] - Vooov_2i_3k[l,c]*T2_1j[l,b,a]  # jik bac

                        V[a,b,c]  = W[a,b,c] + Vovov_1j_3k[b,c]*T1_1i[a] + Vovov_1i_3k[a,c]*T1_1j[b] + Vovov_1i_3j[a,b]*T1_1k[c]
                    end

                    # Compute Energy contribution
                    δjk = j == k
                    for a in 1:v
                        Dijka = Dijk - fv[a]
                        for b in 1:a
                            Dijkab = Dijka - fv[b]
                            δab = a == b
                            for c in 1:b
                                Dd = Dijkab - fv[c]
                                δbc = b == c
                                @inbounds X = (W[a,b,c]*V[a,b,c] + W[a,c,b]*V[a,c,b] + W[b,a,c]*V[b,a,c] + W[b,c,a]*V[b,c,a] + W[c,a,b]*V[c,a,b] + W[c,b,a]*V[c,b,a])
                                @inbounds Y = (V[a,b,c] + V[b,c,a] + V[c,a,b])
                                @inbounds Z = (V[a,c,b] + V[b,a,c] + V[c,b,a])
                                @inbounds Ef = (Y - 2*Z)*(W[a,b,c] + W[b,c,a] + W[c,a,b]) + (Z - 2*Y)*(W[a,c,b]+W[b,a,c]+W[c,b,a]) + 3*X
                                Evals[Threads.threadid()] += Ef*(2-δij-δjk)/(Dd*(1+δab+δbc))
                            end
                        end
                    end
                end
            end
        end
    end

    Et = sum(Evals)
    output("Final (T) contribution: {:15.10f}", Et)
    output("CCSD(T) energy:         {:15.10f}", Et+ccsd.energy)
    return RCCSDpT{T}(ccsd, Et+ccsd.energy, Et)
end
