using LinearAlgebra
using LoopVectorization
using Octavian

function RCCSDpT(Alg::ijk2)
    val = Options.get("return_ints")
    Options.set("return_ints", true)
    ccsd, moints = RCCSD()
    Options.set("return_ints", val)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(mol::Molecule, Alg::ijk2)
    val = Options.get("return_ints")
    Options.set("return_ints", true)
    ccsd, moints = RCCSD(mol)
    Options.set("return_ints", val)
    return RCCSDpT(ccsd, moints, Alg)
end

function RCCSDpT(ccsd::RCCSD, moints::IntegralHelper{T,E,O}, Alg::ijk2) where {T<:AbstractFloat, E<:AbstractERI, O<:AbstractRestrictedOrbitals}

    output("\n   • Perturbative Triples Started\n")
    output("   - Contraction Engine: Octavian")

    T1 = permutedims(ccsd.T1, (2,1))
    T2 = permutedims(ccsd.T2, (4,3,2,1))
    v,o = size(T1)

    # Switch to physicist notation for better memory layout
    ovoo = permutedims(moints["OOOV"], (1,4,2,3))
    vvoo = permutedims(moints["OVOV"], (2,4,1,3))
    vvvo = reshape(permutedims(moints["OVVV"], (3,2,4,1)), v*v, v, o)

    fo = moints["Fii"]
    fv = moints["Faa"]

    # Pre-allocate Intermediate arrays
    Ws  = [Array{T}(undef, v,v,v) for _ = 1:Threads.nthreads()] 
    Vs  = [Array{T}(undef, v,v,v) for _ = 1:Threads.nthreads()] 
    tVs  = [Array{T}(undef, v,v,v) for _ = 1:Threads.nthreads()] 
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
        tV = tVs[id]

        # T1 views → V
        @views T1_i   = T1[:,i]

        # T2 views → W
        @views T2_i   = reshape(T2[:,:,:,i], v*v, o)

        @views vvvo_i = vvvo[:,:,i]

        for j in 1:i
            Dij = fo[i] + fo[j]

            @views vvvo_j = vvvo[:,:,j]

            # T1 views → V
            @views T1_j  = T1[:,j]

            # T2 views → W
            @views T2_j   = reshape(T2[:,:,:,j], v*v, o)
            ## Twice 
            @views T2_ji = T2[:,:,j,i]
            @views T2_ij = T2[:,:,i,j]


            @views vvoo_ij = vvoo[:,:,i,j]

            @views ovoo_ji = ovoo[:,:,j,i]
            @views ovoo_ij = ovoo[:,:,i,j]

            δij = i == j
            for k in 1:j
                Dijk = Dij + fo[k]

                # T1 views → V
                @views T1_k = T1[:,k]
                @views vvvo_k = vvvo[:,:,k]

                # T2 views → W
                @views T2_k   = reshape(T2[:,:,:,k], v*v, o)
                ## Twice 
                @views T2_jk = T2[:,:,j,k]
                @views T2_ik = T2[:,:,i,k]
                @views T2_ki = T2[:,:,k,i]
                @views T2_kj = T2[:,:,k,j]

                @views vvoo_jk = vvoo[:,:,j,k]
                @views vvoo_ik = vvoo[:,:,i,k]
                @views ovoo_ik = ovoo[:,:,i,k]
                @views ovoo_jk = ovoo[:,:,j,k]
                @views ovoo_ki = ovoo[:,:,k,i]
                @views ovoo_kj = ovoo[:,:,k,j]

                # Build intermediates W and V for the set i j k

                # abc
                W = reshape(W, v^2, v)
                matmul_serial!(W, vvvo_j, T2_ik)
                matmul_serial!(W, T2_j, ovoo_ik, -1.0, 1.0)
                W = reshape(W, v,v,v)

                # acb 
                V = reshape(V, v^2, v)
                matmul_serial!(V, vvvo_k, T2_ij)
                matmul_serial!(V, T2_k, ovoo_ij, -1.0, 1.0)
                V = reshape(V, v,v,v)
                permutedims!(tV, V, (1,3,2))
                W .+= tV

                # bac
                V = reshape(V, v^2, v)
                matmul_serial!(V, vvvo_i, T2_jk)
                matmul_serial!(V, T2_i, ovoo_jk, -1.0, 1.0)
                V = reshape(V, v,v,v)
                permutedims!(tV, V, (2,1,3))
                W .+= tV

                # bca
                V = reshape(V, v^2, v)
                matmul_serial!(V, vvvo_k, T2_ji)
                matmul_serial!(V, T2_k, ovoo_ji, -1.0, 1.0)
                V = reshape(V, v,v,v)
                permutedims!(tV, V, (3,1,2))
                W .+= tV

                # cab
                V = reshape(V, v^2, v)
                matmul_serial!(V, vvvo_i, T2_kj)
                matmul_serial!(V, T2_i, ovoo_kj, -1.0, 1.0)
                V = reshape(V, v,v,v)
                permutedims!(tV, V, (2,3,1))
                W .+= tV

                # cba
                V = reshape(V, v^2, v)
                matmul_serial!(V, vvvo_j, T2_ki)
                matmul_serial!(V, T2_j, ovoo_ki, -1.0, 1.0)
                V = reshape(V, v,v,v)
                permutedims!(tV, V, (3,2,1))
                W .+= tV

                @tensor begin
                    V[a,b,c]  = W[a,b,c] + T1_i[a]*vvoo_jk[b,c] + vvoo_ik[a,c]*T1_j[b] + vvoo_ij[a,b]*T1_k[c]
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