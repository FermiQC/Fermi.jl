function RCCSDpT{T}(ccsd::RCCSD, moints::IntegralHelper) where T <: AbstractFloat

    @output "\n   • Perturbative Triples Started\n\n"

    T1 = ccsd.T1.data
    T2 = ccsd.T2.data

    Vvvvo = ints["OVVV"]
    Vvooo = ints["OOOV"]
    Vvovo = ints["OOVV"]

    o,v = size(T1)
    Et::T = 0.0

    fo = diag(ints["FOO"])
    fv = diag(ints["FVV"])

    # Pre-allocate Intermediate arrays
    W  = Array{T}(undef, v,v,v)
    V  = Array{T}(undef, v,v,v)

    # Pre-allocate Views
    Vvvvo_4k =    Array{T}(undef,v,v,v)
    Vvooo_2k_3j = Array{T}(undef,v,o)
    Vvooo_2k_3i = Array{T}(undef,v,o)
    Vvooo_2i_3k = Array{T}(undef,v,o)
    Vvooo_2j_3k = Array{T}(undef,v,o)
    T2_1k_2j    = Array{T}(undef,v,v)
    T2_1k_2i    = Array{T}(undef,v,v) 
    T2_1i_2k    = Array{T}(undef,v,v) 
    T2_1j_2k    = Array{T}(undef,v,v) 
    T2_1k       = Array{T}(undef,o,v,v) 
    T1_1k       = Array{T}(undef,v) 
    Vvovo_2j_4k = Array{T}(undef,v,v) 
    Vvovo_2i_4k = Array{T}(undef,v,v) 
    Vvvvo_4i    = Array{T}(undef,v,v,v) 
    T2_1i       = Array{T}(undef,o,v,v) 
    T1_1i       = Array{T}(undef,v) 
    Vvvvo_4j    = Array{T}(undef,v,v,v) 
    Vvooo_2i_3j = Array{T}(undef,v,o) 
    Vvooo_2j_3i = Array{T}(undef,v,o) 
    T2_1i_2j    = Array{T}(undef,v,v) 
    T2_1j_2i    = Array{T}(undef,v,v) 
    T2_1j       = Array{T}(undef,o,v,v) 
    Vvovo_2i_4j = Array{T}(undef,v,v) 
    T1_1j       = Array{T}(undef,v) 
    @output "Computing energy contribution from occupied orbitals:\n"
    for i in 1:o
        @output "→ Orbital {} of {}\n" i o
        Vvvvo_4i .= view(Vvvvo, :,:,:,i)
        T2_1i    .= view(T2, i,:,:,:)
        T1_1i    .= view(T1, i, :)
        for j in 1:i
            Vvvvo_4j    .= view(Vvvvo, :,:,:,j)
            Vvooo_2i_3j .= view(Vvooo,:,i,j,:)
            Vvooo_2j_3i .= view(Vvooo,:,j,i,:)
            T2_1i_2j    .= view(T2, i,j,:,:)
            T2_1j_2i    .= view(T2, j,i,:,:)
            T2_1j       .= view(T2, j,:,:,:)
            Vvovo_2i_4j .= view(Vvovo, :,i,:,j)
            T1_1j       .= view(T1, j, :)
            δij = Int(i == j)
            for k in 1:j
                Vvvvo_4k    .= view(Vvvvo, :,:,:,k)
                Vvooo_2k_3j .= view(Vvooo,:,k,j,:)
                Vvooo_2k_3i .= view(Vvooo,:,k,i,:)
                Vvooo_2i_3k .= view(Vvooo,:,i,k,:)
                Vvooo_2j_3k .= view(Vvooo,:,j,k,:)

                T2_1k_2j .= view(T2, k,j,:,:)
                T2_1k_2i .= view(T2, k,i,:,:)
                T2_1i_2k .= view(T2, i,k,:,:)
                T2_1j_2k .= view(T2, j,k,:,:)

                T2_1k .= view(T2, k,:,:,:)
                T1_1k .= view(T1, k, :)

                Vvovo_2j_4k .= view(Vvovo, :,j,:,k)
                Vvovo_2i_4k .= view(Vvovo, :,i,:,k)

                # Build intermediates W and V for the set i j k
                @tensoropt begin
                    W[a,b,c]  = (Vvvvo_4i[b,d,a]*T2_1k_2j[c,d] - Vvooo_2k_3j[c,l]*T2_1i[l,a,b]  # ijk abc
                              +  Vvvvo_4i[c,d,a]*T2_1j_2k[b,d] - Vvooo_2j_3k[b,l]*T2_1i[l,a,c]  # ikj acb
                              +  Vvvvo_4k[a,d,c]*T2_1j_2i[b,d] - Vvooo_2j_3i[b,l]*T2_1k[l,c,a]  # kij cab
                              +  Vvvvo_4k[b,d,c]*T2_1i_2j[a,d] - Vvooo_2i_3j[a,l]*T2_1k[l,c,b]  # kji cba
                              +  Vvvvo_4j[c,d,b]*T2_1i_2k[a,d] - Vvooo_2i_3k[a,l]*T2_1j[l,b,c]  # jki bca
                              +  Vvvvo_4j[a,d,b]*T2_1k_2i[c,d] - Vvooo_2k_3i[c,l]*T2_1j[l,b,a]) # jik bac

                    V[a,b,c]  = W[a,b,c] + Vvovo_2j_4k[b,c]*T1_1i[a] + Vvovo_2i_4k[a,c]*T1_1j[b] + Vvovo_2i_4j[a,b]*T1_1k[c]
                end

                # Compute Energy contribution
                δjk = Int(j == k)
                for a in 1:v
                    for b in 1:a
                        δab = Int(a == b)
                        for c in 1:b
                            Dd = fo[i] + fo[j] + fo[k] - fv[a] - fv[b] - fv[c]
                            δbc = Int(b == c)
                            @inbounds X = (W[a,b,c]*V[a,b,c] + W[a,c,b]*V[a,c,b] + W[b,a,c]*V[b,a,c] + W[b,c,a]*V[b,c,a] + W[c,a,b]*V[c,a,b] + W[c,b,a]*V[c,b,a])
                            @inbounds Y = (V[a,b,c] + V[b,c,a] + V[c,a,b])
                            @inbounds Z = (V[a,c,b] + V[b,a,c] + V[c,b,a])
                            E = (Y - 2*Z)*(W[a,b,c] + W[b,c,a] + W[c,a,b]) + (Z - 2*Y)*(W[a,c,b]+W[b,a,c]+W[c,b,a]) + 3*X
                            Et += E*(2-δij-δjk)/(Dd*(1+δab+δbc))
                        end
                    end
                end
            end
        end
    end
    @output "Final (T) contribution: {:15.10f}\n" Et
    @output "CCSD(T) energy:         {:15.10f}\n" Et+ccsd.CorrelationEnergy
    return RCCSDpT{T}(ccsd, Et)
end
