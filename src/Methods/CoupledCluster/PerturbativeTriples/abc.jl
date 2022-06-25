using LinearAlgebra
using LoopVectorization

function RCCSDpT(Alg::abc)
    aoints = IntegralHelper{Float64}()
    return RCCSDpT(aoints, Alg)
end

function RCCSDpT(rhf::RHF, Alg::abc)
    aoints = IntegralHelper()
    moints = IntegralHelper(orbitals=rhf.orbitals)
    ccsd = RCCSD(moints, aoints)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(mol::Molecule, Alg::abc)
    aoints = IntegralHelper{Float64}(molecule=mol)
    return RCCSDpT(aoints, Alg)
end

function RCCSDpT(aoints::IntegralHelper, Alg::abc)
    rhf = RHF(aoints)
    moints = IntegralHelper(orbitals=rhf.orbitals)
    ccsd = RCCSD(moints, aoints)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(ccsd::RCCSD, moints::IntegralHelper{T,E,O}, Alg::abc) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    output("\n   • Perturbative Triples Started\n")
    output("   - Using abc algorithm")

    T1 = ccsd.T1
    T2 = ccsd.T2 

    # Use physicists notation
    ovvv = permutedims(moints["OVVV"], (1,3,2,4))
    ooov = permutedims(moints["OOOV"], (1,3,2,4))
    oovv = permutedims(moints["OVOV"], (1,3,2,4))

    o,v = size(T1)

    fo = moints["Fii"]
    fv = moints["Faa"]

    # Pre-allocate Intermediate arrays
    Ws  = [Array{T}(undef, o,o,o) for _ = 1:Threads.nthreads()] 
    Vs  = [Array{T}(undef, o,o,o) for _ = 1:Threads.nthreads()] 
    Evals = zeros(T, Threads.nthreads())

    output("Computing energy contribution from occupied orbitals:")
    BLAS_THREADS = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    t = @elapsed begin
    Threads.@threads for a in 1:v
    @inbounds begin
        id = Threads.threadid()
        W = Ws[id]
        V = Vs[id]

        # T1 views → V
        @views T1_a   = T1[:,a]

        # T2 views → W
        @views T2_a   = T2[:,:,:,a]

        # ooov views
        @views ooov_a = ooov[:,:,:,a]

        for b in 1:a
            Dab = fv[a] + fv[b]

            # T1 views → V
            @views T1_b  = T1[:,b]

            # T2 views → W
            @views T2_b  = T2[:,:,:,b]
            @views T2_ab = T2[:,:,a,b]
            @views T2_ba = T2[:,:,b,a]

            # ooov views
            @views ooov_b = ooov[:,:,:,b]

            # oovv views
            @views oovv_ab = oovv[:,:,a,b]
            @views oovv_ba = oovv[:,:,b,a]

            # ovvv views
            @views ovvv_ab = ovvv[:,:,a,b]
            @views ovvv_ba = ovvv[:,:,b,a]

            δab = a == b
            for c in 1:b
                Dabc = Dab + fv[c]

                # T1 views → V
                @views T1_c = T1[:,c]

                # T2 views → W
                @views T2_c = T2[:,:,:,c]
                @views T2_ac = T2[:,:,a,c]
                @views T2_ca = T2[:,:,c,a]
                @views T2_bc = T2[:,:,b,c]
                @views T2_cb = T2[:,:,c,b]

                # ooov views
                @views ooov_c = ooov[:,:,:,c]

                # oovv views
                @views oovv_ac = oovv[:,:,a,c]
                @views oovv_ca = oovv[:,:,c,a]
                @views oovv_bc = oovv[:,:,b,c]
                @views oovv_cb = oovv[:,:,c,b]

                # ovvv views
                @views ovvv_ac = ovvv[:,:,a,c]
                @views ovvv_ca = ovvv[:,:,c,a]
                @views ovvv_bc = ovvv[:,:,b,c]
                @views ovvv_cb = ovvv[:,:,c,b]

                # Build intermediates W and V for the set abc
                @tensor begin
                    W[i,j,k]  = ovvv_ab[i,d]*T2_c[j,k,d] - ooov_c[j,k,l]*T2_ab[i,l] #(ijk)(abc)
                    W[i,j,k] += ovvv_ac[i,d]*T2_b[k,j,d] - ooov_b[k,j,l]*T2_ac[i,l] #(ikj)(acb)

                    W[i,j,k] += ovvv_ba[j,d]*T2_c[i,k,d] - ooov_c[i,k,l]*T2_ba[j,l] #(jik)(bac)
                    W[i,j,k] += ovvv_bc[j,d]*T2_c[k,i,d] - ooov_a[k,i,l]*T2_bc[j,l] #(jki)(bca)

                    W[i,j,k] += ovvv_ca[k,d]*T2_b[i,j,d] - ooov_b[i,j,l]*T2_ca[k,l] #(kij)(cab)
                    W[i,j,k] += ovvv_cb[k,d]*T2_c[j,i,d] - ooov_a[j,i,l]*T2_cb[k,l] #(kji)(cba)

                    V[i,j,k] = W[i,j,k] + oovv_bc[j,k]*T1_a[i] + oovv_ac[i,k]*T1_b[j] + oovv_ab[i,j]*T1_c[k]
                end

                # Compute Energy contribution
                δbc = b === c
                for i in 1:o
                    Diabc = fo[i] - Dabc
                    for j in 1:i
                        Dijabc = fo[j] + Diabc
                        δij = i === j
                        @fastmath for k in 1:j
                            Dd = Dijabc + fo[k]
                            δjk = j === k
                            X = (W[i,j,k]*V[i,j,k] + W[i,k,j]*V[i,k,j] + W[j,i,k]*V[j,i,k] + W[j,k,i]*V[j,k,i] + W[k,i,j]*V[k,i,j] + W[k,j,i]*V[k,j,i])
                            Y = (V[i,j,k] + V[j,k,i] + V[k,i,j])
                            Z = (V[i,k,j] + V[j,i,k] + V[k,j,i])
                            Ef  = (Y - 2*Z)*(W[i,j,k] + W[k,i,j] + W[j,k,i]) 
                            Ef += (Z - 2*Y)*(W[i,k,j] + W[j,i,k] + W[k,j,i]) + 3*X
                            Evals[id] += 2*Ef / (Dd * (1+ δij + δjk)*(1 + δab + δbc))
                        end
                    end
                end 
            end
        end
        #output("  Orbital {} ✔️", a)
    end # Threads (i loop)
    end # inbounds
    end # time

    BLAS.set_num_threads(BLAS_THREADS)
    Et = sum(Evals)
    output("Finished in {:5.5f} s", t)
    output("Final (T) contribution: {:15.10f}", Et)
    output("CCSD(T) energy:         {:15.10f}", Et+ccsd.energy)
    return RCCSDpT{T}(ccsd, Et+ccsd.energy, Et)
end