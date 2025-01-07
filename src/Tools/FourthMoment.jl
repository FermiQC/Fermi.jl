using Fermi.Orbitals
using TensorOperations
using ForwardDiff
using ExponentialUtilities


function myfuckingexp(A)
    old = deepcopy(A)
    dA = 1.0

    out = similar(A)
    #out .= [1.0 if i==j]

    fac = 1.0
    Ap = deepcopy(A)
    k = 1.0
    while dA > 1e-6
        old .= out

        fac = fac/k
        Ap *= A
        out += fac * Ap

        println("fac $fac k $k")
        display(out)

        dA = √(sum((out .- old).^2))
        k += 1
    end

    return out
end

function jexp(A)
    outA = deepcopy(A)
    return exponential!(outA, ExponentialUtilities.ExpMethodGeneric())
end

function apply_rotation(κ, C, No, Nv, p =false)

    # Build A from κ
    A = similar(κ, size(C))
    A .= 0.0

    k = 1
    for i = 2:No
        for j = 1:(i-1)
            A[i,j] = κ[k]
            k += 1
        end
    end

    for a = (No+2):(Nv+No)
        for b = (No+1):(a-1)
            A[a,b] = κ[k]
            k += 1
        end
    end

    A = A - A'
    p ? display(jexp(A)) : nothing

    return C*jexp(A)
end

function fourth_moment_localize(orbs)

    bset = BasisSet(orbs.basis, orbs.molecule.atoms)
    println("Getting Integrals")
    t = @elapsed begin
    d4 = GaussianBasis.hexadecapole(bset)
    d3 = GaussianBasis.octupole(bset)
    d2 = GaussianBasis.quadrupole(bset)
    d1 = GaussianBasis.dipole(bset)
    end
    println("Time: $t")

    C = Array(orbs.C)

    No = orbs.molecule.Nα
    Nv = size(C,1) - No

    # Initial Guess
    κ = zeros((No^2 - No) ÷ 2 + (Nv^2 - Nv) ÷ 2)
    println(length(κ))
    println("Initial: $(ξ4(apply_rotation(κ, C, No, Nv), d1, d2, d3, d4, 2))")

    newC = 1
    for _ = 1:10
        t = @elapsed g = get_grad(κ, C, d1, d2, d3, d4, No, Nv)
        println("Time for gradient: $t")
        t = @elapsed h = get_hessian(κ, C, d1, d2, d3, d4, No, Nv)
        println("Time for Hessian: $t")
        κ -= (h^-1) * g
        newC = apply_rotation(κ, C, No, Nv, false)
        println("Final: $(ξ4(newC, d1, d2, d3, d4, 2))")
    end
    return Fermi.Orbitals.GeneralRestrictedOrbitals(orbs.molecule, orbs.basis, 0.0, newC)
end

function get_grad(κ, C, d1, d2, d3, d4, No, Nv)

    #f = A -> ξ4(jexp(A-A')*C, d1, d2, d3, d4, 2)
    f = κ -> ξ4(apply_rotation(κ, C, No, Nv), d1, d2, d3, d4, 2)

    g = x -> ForwardDiff.gradient(f, x)

    return g(κ)
end

function get_hessian(κ, C, d1, d2, d3, d4, No, Nv)

    #f = A -> ξ4(jexp(A-A')*C, d1, d2, d3, d4, 2)
    f = κ -> ξ4(apply_rotation(κ,C, No, Nv), d1, d2, d3, d4, 2)

    h = x -> ForwardDiff.hessian(f, x)

    return h(κ)
end

function ξ4(C, d1, d2, d3, d4, m)

    #@tensoropt begin
    #    dmo1[i,j,x]       := C[μ,i]*C[ν,j]*d1[μ,ν,x]      
    #    dmo2[i,j,x,y]     := C[μ,i]*C[ν,j]*d2[μ,ν,x,y]    
    #    dmo3[i,j,x,y,z]   := C[μ,i]*C[ν,j]*d3[μ,ν,x,y,z]  
    #    dmo4[i,j,x,y,z,w] := C[μ,i]*C[ν,j]*d4[μ,ν,x,y,z,w]
    #end
    #return sum(μ4(dmo1, dmo2, dmo3, dmo4).^m)

    Nao = size(d1, 1)

    Norb    = size(C,2)
    d1_x    = similar(C, Norb, 3)
    d2_xx   = similar(C, Norb, 3)
    d3_xxx  = similar(C, Norb, 3)
    d4_xxxx = similar(C, Norb, 3)

    d2_xy   = similar(C, Norb, 3, 3)
    d3_xxy  = similar(C, Norb, 3, 3)
    d4_xxyy = similar(C, Norb, 3, 3)

    d1_x    .= 0.0
    d2_xx   .= 0.0
    d3_xxx  .= 0.0
    d4_xxxx .= 0.0

    d2_xy   .= 0.0
    d3_xxy  .= 0.0
    d4_xxyy .= 0.0    

    for μ = 1:Nao
        for ν = 1:Nao
            for p = 1:Norb
                Cμν = C[μ,p]*C[ν,p]
                for x = 1:3
                    d1_x[p,x]     += Cμν*d1[μ,ν,x]      
                    d2_xx[p,x]    += Cμν*d2[μ,ν,x,x]    
                    d3_xxx[p,x]   += Cμν*d3[μ,ν,x,x,x]  
                    d4_xxxx[p,x]  += Cμν*d4[μ,ν,x,x,x,x]

                    for y = 1:3
                        d2_xy[p,x,y]    += Cμν*d2[μ,ν,x,y]    
                        d3_xxy[p,x,y]   += Cμν*d3[μ,ν,x,x,y]  
                        d4_xxyy[p,x,y]  += Cμν*d4[μ,ν,x,x,y,y]
                    end
                end
            end
        end
    end

    stuff = μ4(d1_x, d2_xx, d2_xy, d3_xxx, d3_xxy, d4_xxxx, d4_xxyy, 1)
    return sum(stuff.^m)
    #return sum(μ4(C, d1, d2, d3, d4).^m)
end

function μ4(d1, d2, d3, d4)

    norbs = size(d1, 1)
    out = similar(d1, (norbs))
    out .= 0

    for p = eachindex(out)
        for i = 1:3
            out[p] +=  d4[p,p,i,i,i,i]             #   ⟨p|x⁴|p⟩
            out[p] += -4.0*d3[p,p,i,i,i]*d1[p,p,i] # -4⟨p|x³|p⟩⟨p|x|p⟩
            out[p] +=  6.0*d2[p,p,i,i]*d1[p,p,i]^2 #  6⟨p|x²|p⟩⟨p|x|p⟩²
            out[p] += -3.0*d1[p,p,i]^4             # -3⟨p|x|p⟩⁴
            out[p] += -3.0*d1[p,p,i]^4             # -3⟨p|x|p⟩⁴
            for j = 1:(i-1) # i > j
                out[p] +=  2.0*d4[p,p,i,i,j,j]                 #  2⟨p|xᵢ²xⱼ²|p⟩
                out[p] += -4.0*d3[p,p,i,i,j]*d1[p,p,j]         # -4⟨p|xᵢ²xⱼ|p⟩⟨p|xⱼ|p⟩
                out[p] += -4.0*d3[p,p,j,j,i]*d1[p,p,i]         # -4⟨p|xⱼ²xᵢ|p⟩⟨p|xᵢ|p⟩
                out[p] +=  2.0*d2[p,p,i,i]*d1[p,p,j]^2         #  2⟨p|xᵢ²|p⟩⟨p|xⱼ|p⟩²
                out[p] +=  2.0*d2[p,p,j,j]*d1[p,p,i]^2         #  2⟨p|xⱼ²|p⟩⟨p|xᵢ|p⟩²
                out[p] += -6.0*(d1[p,p,j]*d1[p,p,i])^2         # -6⟨p|xⱼ|p⟩²⟨p|xᵢ|p⟩²
                out[p] +=  8.0*d2[p,p,i,j]*d1[p,p,i]*d1[p,p,j] #  8⟨p|xᵢxⱼ|p⟩⟨p|xᵢ|p⟩⟨p|xⱼ|p⟩
            end
        end
    end

    return out
end

function μ4(d_x, d_xx, d_xy, d_xxx, d_xxy, d_xxxx, d_xxyy, x::Int)

    norbs = size(d_x, 1)
    out = similar(d_x, (norbs))
    out .= 0

    for p = eachindex(out)
        for i = 1:3
            out[p] +=  d_xxxx[p,i]                #   ⟨p|x⁴|p⟩
            out[p] += -4.0*d_xxx[p,i]*d_x[p,i]    # -4⟨p|x³|p⟩⟨p|x|p⟩
            out[p] +=  6.0*d_xx[p,i]*d_x[p,i]^2   #  6⟨p|x²|p⟩⟨p|x|p⟩²
            out[p] += -3.0*d_x[p,i]^4             # -3⟨p|x|p⟩⁴
            for j = 1:(i-1) # i > j
                out[p] +=  2.0*d_xxyy[p,i,j]                 #  2⟨p|xᵢ²xⱼ²|p⟩
                out[p] += -4.0*d_xxy[p,i,j]*d_x[p,j]        # -4⟨p|xᵢ²xⱼ|p⟩⟨p|xⱼ|p⟩
                out[p] += -4.0*d_xxy[p,j,i]*d_x[p,i]        # -4⟨p|xⱼ²xᵢ|p⟩⟨p|xᵢ|p⟩
                out[p] +=  2.0*d_xx[p,i]*d_x[p,j]^2          #  2⟨p|xᵢ²|p⟩⟨p|xⱼ|p⟩²
                out[p] +=  2.0*d_xx[p,j]*d_x[p,i]^2        #  2⟨p|xⱼ²|p⟩⟨p|xᵢ|p⟩²
                out[p] += -6.0*(d_x[p,j]*d_x[p,i])^2         # -6⟨p|xⱼ|p⟩²⟨p|xᵢ|p⟩²
                out[p] +=  8.0*d_xy[p,i,j]*d_x[p,i]*d_x[p,j] #  8⟨p|xᵢxⱼ|p⟩⟨p|xᵢ|p⟩⟨p|xⱼ|p⟩
            end
        end
    end

    return out
end