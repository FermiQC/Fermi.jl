using LinearAlgebra

function RCCSDpT(rhf::RHF, Alg::RpTa)
    aoints = IntegralHelper()
    moints = IntegralHelper(orbitals=rhf.orbitals)
    ccsd = RCCSD(moints, aoints)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(Alg::RpTa)
    aoints = IntegralHelper{Float64}()
    rhf = RHF(aoints)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    ccsd = RCCSD(moints, aoints)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(ccsd::RCCSD, moints::IntegralHelper{T,E,O}, Alg::RpTa) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    output("\n   • Perturbative Triples Started\n")

    T1 = ccsd.T1.data
    T2 = ccsd.T2.data

    Vovvv = moints["OVVV"].data
    Vooov = moints["OOOV"].data
    Vovov = moints["OVOV"].data

    o,v = size(T1)

    fo = moints["Fii"].data
    fv = moints["Faa"].data

    # Pre-allocate Intermediate arrays
    Ws  = [Array{T}(undef, v,v,v) for i = 1:Threads.nthreads()] 
    Vs  = [Array{T}(undef, v,v,v) for i = 1:Threads.nthreads()] 
    #Evals = zeros(T, Threads.nthreads())
    Evals = zeros(T, o)

    output("Computing energy contribution from occupied orbitals:")
    BLAS_THREADS = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)
    t = @elapsed begin
    Threads.@threads for i in 1:o
        id = Threads.threadid()
        W = Ws[id]
        V = Vs[id]
        @views Vovvv_1i = Vovvv[i,:,:,:]

        @views T2_1i    = T2[i,:,:,:]
        @views T1_1i    = T1[i, :]
        for j in 1:i
            Dij = fo[i] + fo[j]

            @views Vovvv_1j    = Vovvv[j,:,:,:]
            @views Vovov_1i_3j = Vovov[i,:,j,:]

            @views Vooov_2i_3j = Vooov[:,i,j,:]
            @views Vooov_2j_3i = Vooov[:,j,i,:]

            @views T2_1i_2j    = T2[i,j,:,:]
            @views T2_1j_2i    = T2[j,i,:,:]
            @views T2_1j       = T2[j,:,:,:]
            @views T1_1j       = T1[j, :]
            δij = i == j
            for k in 1:j
                Dijk = Dij + fo[k]

                @views Vovvv_1k    = Vovvv[k,:,:,:]

                @views Vooov_2j_3k = Vooov[:,j,k,:]
                @views Vooov_2k_3j = Vooov[:,k,j,:]
                @views Vooov_2k_3i = Vooov[:,k,i,:]
                @views Vooov_2i_3k = Vooov[:,i,k,:]

                @views T2_1k_2j = T2[k,j,:,:]
                @views T2_1k_2i = T2[k,i,:,:]
                @views T2_1i_2k = T2[i,k,:,:]
                @views T2_1j_2k = T2[j,k,:,:]

                @views T2_1k = T2[k,:,:,:]
                @views T1_1k = T1[k, :]

                @views Vovov_1j_3k = Vovov[j,:,k,:]
                @views Vovov_1i_3k = Vovov[i,:,k,:]

                # Build intermediates W and V for the set i j k
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
                            @inbounds begin
                                X = (W[a,b,c]*V[a,b,c] + W[a,c,b]*V[a,c,b] + W[b,a,c]*V[b,a,c] + W[b,c,a]*V[b,c,a] + W[c,a,b]*V[c,a,b] + W[c,b,a]*V[c,b,a])
                                Y = (V[a,b,c] + V[b,c,a] + V[c,a,b])
                                Z = (V[a,c,b] + V[b,a,c] + V[c,b,a])
                                Ef = (Y - 2*Z)*(W[a,b,c] + W[b,c,a] + W[c,a,b]) + (Z - 2*Y)*(W[a,c,b]+W[b,a,c]+W[c,b,a]) + 3*X
                            end
                            Evals[i] += Ef*(2-δij-δjk)/(Dd*(1+δab+δbc))
                        end
                    end
                end
            end
        end
        output("  Orbital {} ✔️", i)
    end
    end

    LinearAlgebra.BLAS.set_num_threads(BLAS_THREADS)
    Et = sum(Evals)
    output("Finished in {:5.5f} s", t)
    output("Final (T) contribution: {:15.10f}", Et)
    output("CCSD(T) energy:         {:15.10f}", Et+ccsd.energy)
    return RCCSDpT{T}(ccsd, Et+ccsd.energy, Et)
end
