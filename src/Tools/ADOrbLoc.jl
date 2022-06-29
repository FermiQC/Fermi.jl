using ExponentialUtilities
using TensorOperations
using ForwardDiff

function optimize(C, eri, ndocc, γ = 1.0)
    #m = metric2(C, eri, ndocc)
    m = self_repulsion_metric(C, eri)
    e = exact_metric(C, eri, ndocc)
    println("Initial C-Metric: $m")
    println("Initial E-Metric: $e")

    A = zeros(size(C))

    ite = 1
    while true
        if ite > 10
            break
        end
        println("Iter $ite")
        dA = get_grad(A, C, eri, ndocc)
        #dA = get_exact_grad(A, C, eri, ndocc)

        # Delete OV mixing
        A[1:ndocc, (ndocc+1):end] .= 0.0
        A[(ndocc+1):end, 1:ndocc] .= 0.0

        # Normalize
        maxA = maximum(abs.(A))
        A ./ maxA

        # Apply gradient - Newton method
        for i = 1:ndocc
            for j = 1:ndocc

                #if dA[i,j] < 1e-8
                #    continue
                #end
                A[i,j] = A[i,j] - γ*dA[i,j]
            end
        end

        for a = (ndocc+1):size(C,1) 
            for b = (ndocc+1):size(C,1) 
                #if dA[a,b] < 1e-8
                #    continue
                #end
                A[a,b] = A[a,b] - γ*dA[a,b]
            end
        end

        newC = jexp(A-A')*C 
        #m1 = metric2(newC, eri, ndocc)
        m1 = self_repulsion_metric(newC, eri)
        e = exact_metric(newC, eri, ndocc)
        println("New Metric: $m1")
        println("New Count: $e")
        ite += 1
    end
end

function optimize2(metric, C, eri, ndocc, γ = 1.0; maxite=10)
    m = metric(C, eri, ndocc)
    e = exact_metric(C, eri, ndocc)

    println("2 -Initial C-Metric: $m")
    println("2 - Initial E-Metric: $e")

    ite = 1
    while true
        if ite > maxite
            break
        end
        A = zeros(size(C))
        println("Iter $ite")
        dA = get_grad(metric, A, C, eri, ndocc)

        # Delete OV mixing
        A[1:ndocc, (ndocc+1):end] .= 0.0
        A[(ndocc+1):end, 1:ndocc] .= 0.0

        # Normalize
        #maxA = maximum(abs.(A))
        #A ./ maxA

        # Apply gradient - Newton method
        for i = 1:ndocc
            for j = 1:ndocc

                #if dA[i,j] < 1e-8
                #    continue
                #end
                A[i,j] = A[i,j] - γ*dA[i,j]
            end
        end

        for a = (ndocc+1):size(C,1) 
            for b = (ndocc+1):size(C,1) 
                #if dA[a,b] < 1e-8
                #    continue
                #end
                A[a,b] = A[a,b] - γ*dA[a,b]
            end
        end

        C = jexp(A-A')*C 
        m1 = metric(C, eri, ndocc)
        e = exact_metric(C, eri, ndocc)
        println("New Metric: $m1")
        println("New Count: $e")
        ite += 1
    end
end

function get_grad(metric, A0, C, eri, ndocc)

    f = A -> metric(jexp(A-A')*C, eri, ndocc)

    g = x -> ForwardDiff.gradient(f, x)

    return g(A0)
end

function get_exact_grad(A0, C, eri, ndocc)

    f = A -> exact_metric(jexp(A-A')*C, eri, ndocc)

    return grad(central_fdm(5,1), f, A0)[1]
end

function metric2(C, eri, ndocc)
    o = 1:ndocc
    v = (ndocc+1):size(C,1)

    @views Co = C[:,o]
    @views Cv = C[:,v]

    @tensoropt iajb[i,a,j,b] := Co[μ,i]*Cv[ν,a]*Co[λ,j]*Cv[σ,b]*eri[μ,ν,λ,σ]

    #av = sum(iajb) / length(iajb)
    #S2 = sum((iajb .- av).^2) / (length(iajb) - 1)
    #return S2
    mv = maximum(abs.(iajb))
    return sum(act.(iajb ./ mv, 1)) 
    #return sum(iajb.^2) / length(iajb)
    #return sum(log.(1.0 .+ iajb.^2)) / length(iajb)
end

function exact_metric(C, eri, ndocc)
    o = 1:ndocc
    v = (ndocc+1):size(C,1)

    @views Co = C[:,o]
    @views Cv = C[:,v]

    @tensoropt iajb[i,a,j,b] := Co[μ,i]*Cv[ν,a]*Co[λ,j]*Cv[σ,b]*eri[μ,ν,λ,σ]

    return count(x-> abs(x) > 1e-8, iajb) / length(iajb)
end

function jexp(A)
    outA = deepcopy(A)
    return exponential!(outA, ExponentialUtilities.ExpMethodGeneric())
end

function act(x, ϵ = 100)
    return tanh(ϵ*x)^2
end

function self_repulsion_metric(C, eri, ndocc)
    out = 0.0
    for p = 1:size(C,1)
        @views Cp = C[:,p]
        @tensoropt I = Cp[μ]*Cp[ν]*Cp[ρ]*Cp[σ]*eri[μ,ν,ρ,σ]
        out += I
    end
    return out
end

function Cdensity(C, eri, ndocc)
   return sum(act.(C, 10)) / length(C)
end