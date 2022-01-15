using Fermi.Orbitals
using TensorOperations

function boys_localize(orbs::AbstractRestrictedOrbitals)
    nothing
end

function ER_localize(orbs::AbstractRestrictedOrbitals, ints, mix_range)

    eri = ints["ERI"]

    # Maximaze the self-repulsion (pp|pp)
    newC = deepcopy(orbs.C)

    N = length(mix_range)

    println("Number of mixing pairs: $N")
    # Create a list of pairs to be rotated
    Npairs = N*(N-1) ÷ 2
    Δvals = zeros(Npairs)
    θvals = zeros(Npairs)
    ij_vals = Array{Tuple{Int64, Int64}}(undef, Npairs)

    ## Fill ij_vals
    k = 1
    for i in mix_range
        for j in mix_range
            if i > j
                ij_vals[k] = (i,j)
                k += 1
            end
        end
    end

    #@views Co = newC[:, o]
    #self = zeros(No)
    #J = zeros(No, No)
    #K = zeros(No, No)
    #L = zeros(No, No)
    #for p = o
    #    @views C1 = newC[:,p]
    #    @tensoropt A = C1[μ]*C1[ν]*C1[λ]*C1[σ]*eri[μ, ν, λ, σ]
    #    self[p] = A
    #    for q = o
    #        @views C2 = newC[:,q]
    #        @tensoropt B = C1[μ]*C1[ν]*C2[λ]*C2[σ]*eri[μ, ν, λ, σ]
    #        @tensoropt C = C1[μ]*C2[ν]*C1[λ]*C2[σ]*eri[μ, ν, λ, σ]
    #        @tensoropt D = C1[μ]*C1[ν]*C1[λ]*C2[σ]*eri[μ, ν, λ, σ]

    #        J[p,q] = B
    #        K[p,q] = C
    #        L[p,q] = D
    #    end
    #end

    println("Self repulsion: $(selfrepulsion(newC, mix_range, ints))")

    ite = 1
    for idx in eachindex(ij_vals)
        (i,j) = ij_vals[idx]

        #iiii = self[i]
        #jjjj = self[j]

        #iiij = L[i,j]
        #jjij = L[j,i]

        #iijj = J[i,j]
        #ijij = K[i,j]

        #B = iiij - jjij
        #A = ijij + 0.5*iijj - 0.25*(iiii + jjjj)
        #θ = 0.25 * atan(-B / A)

        #Δvals[idx] = A*(1 - cos(4θ)) + B*sin(4θ)
        #θvals[idx] = θ

        θvals[idx], Δvals[idx] = ER_rotation_angle(newC, i,j, ints["ERI"])
    end

    while true
        if ite > 30
            println("Ops did not converge")
            break
        end

        println("Iter $ite")

        RMS = √(sum((Δvals) .^2) / length(Δvals))

        println("ΔI RMS: $RMS")

        if RMS < 1e-5
            break
        end
        
        _, idx = findmax(Δvals)

        i,j = ij_vals[idx]
        θ = θvals[idx]
        Δ = Δvals[idx]

        if Δ < 1e-8
            println("Converged")
            break
        end

        println("Rotating $i and $j by $θ. Δ expected: $Δ")
        rotate_orbitals!(newC, θ, i,j)
        θvals[idx], Δvals[idx] = ER_rotation_angle(newC, i,j, ints["ERI"])

        #@assert abs(Θvals[idx]) < 1e-10
        #@assert abs(Δvals[idx]) < 1e-10

        println("Updating integrals")

        for k in mix_range
            # Update i vals
            if k > i
                idx, = findall(x->x==(k,i), ij_vals)  # PROVISORIO
                θvals[idx], Δvals[idx] = ER_rotation_angle(newC, k,i, ints["ERI"])
            elseif i > k
                idx, = findall(x->x==(i,k), ij_vals)  # PROVISORIO
                θvals[idx], Δvals[idx] = ER_rotation_angle(newC, i,k, ints["ERI"])
            end

            # Update j vals
            if k > j
                idx, = findall(x->x==(k,j), ij_vals)  # PROVISORIO
                θvals[idx], Δvals[idx] = ER_rotation_angle(newC, k,j, ints["ERI"])
            elseif j > k
                idx, = findall(x->x==(j,k), ij_vals)  # PROVISORIO
                θvals[idx], Δvals[idx] = ER_rotation_angle(newC, j,k, ints["ERI"])
            end
        end

        println("Itotal after rotation: $(selfrepulsion(newC, mix_range, ints))")
        ite += 1
    end

    energy = sd_energy(orbs, orbs.molecule.Nα, ints)
    neworbs = Fermi.Orbitals.GeneralRestrictedOrbitals(orbs.molecule, orbs.basis, energy, newC)
    return neworbs
end

function selfrepulsion(C, r, ints)
    sum([get_pppp(C, ints["ERI"], i) for i = r])
end

function ER_rotation_angle(C, i,j, eri)

    @views Ci = C[:,i]
    @views Cj = C[:,j]

    iiii = @tensoropt Ci[μ]*Ci[ν]*Ci[ρ]*Ci[σ]*eri[μ,ν,ρ,σ]
    jjjj = @tensoropt Cj[μ]*Cj[ν]*Cj[ρ]*Cj[σ]*eri[μ,ν,ρ,σ]

    iiij = @tensoropt Ci[μ]*Ci[ν]*Ci[λ]*Cj[σ]*eri[μ,ν,λ,σ]
    jjij = @tensoropt Cj[μ]*Cj[ν]*Ci[λ]*Cj[σ]*eri[μ,ν,λ,σ]

    iijj = @tensoropt Ci[μ]*Ci[ν]*Cj[λ]*Cj[σ]*eri[μ,ν,λ,σ]
    ijij = @tensoropt Ci[μ]*Cj[ν]*Ci[λ]*Cj[σ]*eri[μ,ν,λ,σ]

    B = iiij - jjij
    A = ijij + 0.5*iijj - 0.25*(iiii + jjjj)
    θ = 0.25 * atan(-B / A)

    Δ = A*(1 - cos(4θ)) + B*sin(4θ)

    return θ, Δ
end

function rotate_orbitals!(C, θ, i, j)
    ui = C[:,i]
    uj = C[:,j]

    C[:,i] =  ui .* cos(θ) + uj .* sin(θ)
    C[:,j] = -ui .* sin(θ) + uj .* cos(θ)
end

# Get mo (pp|pp) integral
function get_pppp(C, aoeri, p)
    @views Cp = C[:,p]
    @tensoropt I = Cp[μ]*Cp[ν]*Cp[ρ]*Cp[σ]*aoeri[μ,ν,ρ,σ]
    return I
end

# Get mo (pp|pq) integral
function get_pppq(C, aoeri, p, q)
    @views Cp = C[:,p]
    @views Cq = C[:,q]
    @tensoropt I = Cp[μ]*Cp[ν]*Cp[ρ]*Cq[σ]*aoeri[μ,ν,ρ,σ]
    return I
end

# Get mo (pp|qq) integral
function get_ppqq(C, aoeri, p, q)
    @views Cp = C[:,p]
    @views Cq = C[:,q]
    @tensoropt I = Cp[μ]*Cp[ν]*Cq[ρ]*Cq[σ]*aoeri[μ,ν,ρ,σ]
    return I
end

# Get mo (pq|pq) integral
function get_pqpq(C, aoeri, p, q)
    @views Cp = C[:,p]
    @views Cq = C[:,q]
    @tensoropt I = Cp[μ]*Cq[ν]*Cp[ρ]*Cq[σ]*aoeri[μ,ν,ρ,σ]
    return I
end

function sd_energy(orbs, ndocc, ints)

    C = orbs.C
    Co = C[:,1:ndocc]
    D = Co * Co'

    H = ints["T"] + ints["V"]
    F = similar(H)

    Fermi.HartreeFock.build_fock!(F, H, D, ints)

    E = Fermi.HartreeFock.RHFEnergy(D, H, F)

    return E + orbs.molecule.Vnuc
end