using LinearAlgebra
using LoopVectorization

function RCCSDpT(Alg::ijk)
    val = Options.get("return_ints")
    Options.set("return_ints", true)
    ccsd, moints = RCCSD()
    Options.set("return_ints", val)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(mol::Molecule, Alg::ijk)
    val = Options.get("return_ints")
    Options.set("return_ints", true)
    ccsd, moints = RCCSD(mol)
    Options.set("return_ints", val)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(ccsd::RCCSD, moints::IntegralHelper{T,E,O}, Alg::ijk) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    output("\n   • Perturbative Triples Started\n")

    T1 = permutedims(ccsd.T1, (2,1))
    T2 = permutedims(ccsd.T2, (4,3,2,1))

    # Invert order, but retains chemist's notation
    Vvvvo = permutedims(moints["OVVV"], (4,3,2,1))

    # Switch to physicist notation for better memory layout
    Vvooo = permutedims(moints["OOOV"], (4,1,2,3))
    Vvvoo = permutedims(moints["OVOV"], (2,4,1,3))

    v,o = size(T1)

    fo = moints["Fii"]
    fv = moints["Faa"]

    # Pre-allocate Intermediate arrays
    Ws  = [Array{T}(undef, v,v,v) for _ = 1:Threads.nthreads()] 
    Vs  = [Array{T}(undef, v,v,v) for _ = 1:Threads.nthreads()] 
    Evals = zeros(T, Threads.nthreads())

    output("Computing energy contribution from occupied orbitals:")
    BLAS_THREADS = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    t = @elapsed begin
    Threads.@threads for i in 1:o
    @inbounds begin
        id = Threads.threadid()
        W = Ws[id]
        V = Vs[id]

        # T1 views → V
        @views T1_i   = T1[:,i]

        # T2 views → W
        @views T2_i   = T2[:,:,:,i]

        @views Vvvvo_i = Vvvvo[:,:,:,i]

        for j in 1:i
            Dij = fo[i] + fo[j]

            # T1 views → V
            @views T1_j  = T1[:,j]

            # T2 views → W
            @views T2_j  = T2[:,:,:,j]
            ## Twice 
            @views T2_ji = T2[:,:,j,i]
            @views T2_ij = T2[:,:,i,j]

            @views Vvvvo_j    = Vvvvo[:,:,:,j]

            @views Vvvoo_ij = Vvvoo[:,:,i,j]

            @views Vvooo_ij = Vvooo[:,:,i,j]
            @views Vvooo_ji = Vvooo[:,:,j,i]

            δij = i == j
            for k in 1:j
                Dijk = Dij + fo[k]

                # T1 views → V
                @views T1_k = T1[:,k]

                # T2 views → W
                @views T2_k = T2[:,:,:,k]
                ## Twice 
                @views T2_jk = T2[:,:,j,k]
                @views T2_ik = T2[:,:,i,k]
                @views T2_ki = T2[:,:,k,i]
                @views T2_kj = T2[:,:,k,j]

                @views Vvvvo_k  = Vvvvo[:,:,:,k]

                @views Vvooo_jk = Vvooo[:,:,j,k]
                @views Vvooo_kj = Vvooo[:,:,k,j]
                @views Vvooo_ki = Vvooo[:,:,k,i]
                @views Vvooo_ik = Vvooo[:,:,i,k]

                @views Vvvoo_jk = Vvvoo[:,:,j,k]
                @views Vvvoo_ik = Vvvoo[:,:,i,k]

                # Build intermediates W and V for the set i j k
                @tensor begin
                    W[a,b,c] =  Vvvvo_i[d,b,a]*T2_jk[d,c] - Vvooo_ij[b,l]*T2_k[a,c,l]
                    W[a,b,c] += Vvvvo_j[d,a,b]*T2_ik[d,c] - Vvooo_ji[a,l]*T2_k[b,c,l]
                    W[a,b,c] += Vvvvo_k[d,a,c]*T2_ij[d,b] - Vvooo_ki[a,l]*T2_j[c,b,l]
                    W[a,b,c] += Vvvvo_k[d,b,c]*T2_ij[a,d] - Vvooo_kj[b,l]*T2_i[c,a,l]
                    W[a,b,c] += Vvvvo_i[d,c,a]*T2_jk[b,d] - Vvooo_ik[c,l]*T2_j[a,b,l]
                    W[a,b,c] += Vvvvo_j[d,c,b]*T2_ik[a,d] - Vvooo_jk[c,l]*T2_i[b,a,l]

                    V[a,b,c]  = W[a,b,c] + T1_i[a]*Vvvoo_jk[b,c] + Vvvoo_ik[a,c]*T1_j[b] + Vvvoo_ij[a,b]*T1_k[c]
                end

                # Compute Energy contribution
                δjk = j === k
                for a in 1:v
                    Dijka = Dijk - fv[a]
                    for b in 1:a
                        Dijkab = Dijka - fv[b]
                        δab = a === b
                        @fastmath for c in 1:b
                            Dd = Dijkab - fv[c]
                            δbc = b === c
                            X = (W[a,b,c]*V[a,b,c] + W[a,c,b]*V[a,c,b] + W[b,a,c]*V[b,a,c] + W[b,c,a]*V[b,c,a] + W[c,a,b]*V[c,a,b] + W[c,b,a]*V[c,b,a])
                            Y = (V[a,b,c] + V[b,c,a] + V[c,a,b])
                            Z = (V[a,c,b] + V[b,a,c] + V[c,b,a])
                            Ef = (Y - 2*Z)*(W[a,b,c] + W[b,c,a] + W[c,a,b]) + (Z - 2*Y)*(W[a,c,b]+W[b,a,c]+W[c,b,a]) + 3*X
                            Evals[id] += Ef*(2 - δij - δjk) / (Dd*(1 + δab + δbc))
                        end
                    end
                end 
            end
        end
        output("  Orbital {} ✔️", i)
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