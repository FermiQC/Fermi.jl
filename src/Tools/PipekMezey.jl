using Fermi
using LinearAlgebra

function mulliken_charges(Co, ints)

    S = ints["S"]
    q = diag(C0*transpose(C0)*S)

    ## This might change soon
    basis = ints.orbitals.basisset
    natoms = basis.natoms
    basis_per_atom = zeros(Int, natoms)
    charges = zeros(natoms)
    for a in 1:natoms
        n = 0
        for b in basis[a]
            n += 2*b.l + 1
        end
        basis_per_atom[a] = n
    end

    off = 0
    for a in 1:natoms
        n = basis_per_atom[a]
        charges[a] = sum(q[(off+1):(off+n)])
        off += n
    end

    return charges
end

function Qia_mulliken(Co, ints)
    Nocc = size(Co, 2)
    S = ints["S"]

    basis = ints.orbitals.basisset
    natoms = basis.natoms

    basis_per_atom = zeros(Int, natoms)
    for a in 1:natoms
        n = 0
        for b in basis[a]
            n += 2*b.l + 1
        end
        basis_per_atom[a] = n
    end

    Qia = zeros(Nocc, natoms)
    off = 0
    for a in 1:natoms
        n = basis_per_atom[a]
        for i in 1:Nocc
            Ci = Co[:,i]
            q = diag(Ci*transpose(Ci)*S)
            Qia[i,a] += sum(q[(off+1):(off+n)])
        end
        off += n
    end

    return Qia
end

function Qija_mulliken(Co, ints)
    Nocc = size(Co, 2)
    S = ints["S"]

    basis = ints.orbitals.basisset
    natoms = basis.natoms

    basis_per_atom = zeros(Int, natoms)
    for a in 1:natoms
        n = 0
        for b in basis[a]
            n += 2*b.l + 1
        end
        basis_per_atom[a] = n
    end

    Qija = zeros(Nocc, Nocc, natoms)
    off = 0
    for a in 1:natoms
        n = basis_per_atom[a]
        for i in 1:Nocc
            for j in i:Nocc
                Ci = Co[:,i]
                Cj = Co[:,j]
                q = diag(Ci*transpose(Cj)*S + Cj*transpose(Ci)*S)
                Qija[i,j,a] = sum(q[(off+1):(off+n)])
                if i != j
                    Qija[j,i,a] = Qija[i,j,a]
                end
            end
        end
        off += n
    end

    return Qija
end

function pm_delocalization_index(Co, ints)
    Qia2 = Qia_mulliken(Co, ints) .^ 2
    di = sum(Qia2[:,a] for a = 1:size(Qia2,2))
    return 1 ./ di
end

function pm_metric(Co, ints)
    sum(Qia_mulliken(Co, ints) .^ 2)
end

function jacobi_sweep(C, ints)

    Norbs = size(C, 2)
    Natoms = length(ints.molecule.atoms)

    output("Initiating Pipek-Mezey Localization")
    output(" • Algorithm: Jacobi Sweeps (Pair-wise rotations)\n")

    output("Initial PM metric: $(pm_metric(C, ints))")

    θvals = zeros(Norbs, Norbs)

    max_iter = 1000
    iter = 1
    δC = 1.0
    while δC > 1e-10

        if iter > max_iter
            break
        end
        output("\n\n Iter $iter")

        # Get array for Qᵃₖ = ∑(μ∈A) ∑(ν) Cμk⋅Ckν⋅Sνμ   
        Qia = Qia_mulliken(C, ints)

        # Get array for Qᵃᵢⱼ = ∑(μ∈A) ∑(ν) (Cμi⋅Cjν⋅Sνμ + Cμj⋅Ciν⋅Sνμ)
        Qija = Qija_mulliken(C, ints)

        for i = 1:Norbs
            for j = 1:(i-1)
                A = 0.0
                B = 0.0
                for a = 1:Natoms
                    A += 2*Qija[i,j,a]*(Qia[i,a] - Qia[j,a])
                    B += Qija[i,j,a]^2 - (Qia[i,a] - Qia[j,a])^2
                end
                θvals[i,j] = atan(-A/B)/4
                if i != j 
                    θvals[j,i] = θvals[i,j]
                end
            end
        end

        θmax, cmax = findmax(θvals)
        θmin, cmin = findmin(θvals)
        θ = abs(θmax) > abs(θmin) ? θmax : θmin
        i,j = abs(θmax) > abs(θmin) ? cmax.I : cmin.I

        output("Orbitals $i and $j to be rotated by $(rad2deg(θ))")

        # Rotate
        Cnew = copy(C)
        Cnew[:, i] =  cos(θ)*C[:,i] + sin(θ)*C[:,j]
        Cnew[:, j] = -sin(θ)*C[:,i] + cos(θ)*C[:,j]

        output("Final PM metric: $(pm_metric(Cnew, ints))")
        δC = √(sum((Cnew .- C) .^2)/length(C))

        C .= Cnew

        iter += 1
    end
    θvals
end