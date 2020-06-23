module Davidson
using LinearAlgebra
export eigdav
function eigdav(A, eigs, k, kmax, tol)
    n = size(A, 1)
    V = zeros((n, n))
    theta = 0
    w = 0
    theta_old = 0
    I = Diagonal(ones(n, n))
    t = Diagonal(randn(n, n))
    for m = k:k:kmax
        if m <= k
            for j = 1:1:k
                V[:, j] = t[:, j]#/norm(t[:,j])
            end
            theta_old = ones(eigs)
        else
            theta_old = theta[1:eigs]
        end
        F = qr(V)
        V = Array{Float64}(F.Q)
        @views T = transpose(V[:, 1:(m+1)]) * A * V[:, 1:(m+1)]
        theta, s = eigen(T)

        #THETA = eigvals(T)
        #S = eigvecs(T)
        #idx = sortperm(THETA)
        #theta = THETA[idx]
        #s = S[:,idx]
        #println(theta)
        for j = 1:1:k
            @views w = (A - theta[j] * I) * V[:, 1:(m+1)] * s[:, j]
            println("DIFF ", theta[j] - A[j, j])
            q = w / (theta[j] - A[j, j])
            V[:, m + j + 1] = q
        end
        normm = norm(theta[1:eigs] - theta_old)
        println("NORM ", normm)
        if normm < tol
            return theta[1:eigs]
        end

    end
    return theta[1:eigs]
end
end #module Davidson
