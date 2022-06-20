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

    """
    Equations from J. Chem. Phys. 90, 4916 (1989); https://doi.org/10.1063/1.456588
    Metric     -> Equation 18
    Functional -> Equation 28
    Ast        -> Eqaution 29a
    Bst        -> Equation 29b
    cos(4α)    -> Equation 13b
    ⟨s|Pₐ|t⟩   -> Equation 31
    """
    Norbs = size(C, 2)
    Natoms = length(ints.molecule.atoms)

    output("Initiating Pipek-Mezey Localization")
    output(" • Algorithm: Jacobi Sweeps (Pair-wise rotations)\n")

    if Norbs < 2
        output("Number of orbitals must be greater than 1 for localization.")
        return C
    end

    output("Initial PM metric: $(pm_metric(C, ints))")

    θvals = zeros(Norbs, Norbs)

    max_iter = 100
    iter = 0
    δC = 1.0
    Cnew = deepcopy(C)
    Cold = deepcopy(C)

    output("\n\n{:5s}     {:<10s}     {:10s}     {:10s}", "Iter", "Metric", "RMS", "# Rotations")
    metric = pm_metric(Cnew, ints)
    output("{:5d}     {:<10.5f}", iter, metric)
    iter +=1

    while δC > 1e-10

        if iter > max_iter
            break
        end

        # Get ⟨s|Pₐ|t⟩
        P = get_Pst(Cold, ints)

        for s = 1:Norbs
            for t = 1:(s-1)
                Ast = 0.0
                Bst = 0.0
                for a = 1:Natoms
                    Ast += P[s,t,a]^2 - 0.25*(P[s,s,a] - P[t,t,a])^2
                    Bst += P[s,t,a]*(P[s,s,a] - P[t,t,a])
                end
                AB = Ast^2 + Bst^2
                if abs(AB) > 1e-16
                    cos4θ = -Ast / √AB
                    θvals[s,t] = 0.25 * acos(cos4θ) * (Bst > 0 ? 1 : -1)
                else
                    θvals[s,t] = 0.0
                end
            end
        end

        # Rotate angles until all orbitals have been roated
        _, cmax = findmax(abs, θvals)
        θij = θvals[cmax]
        i,j = cmax.I
        nrot = 0
        while abs(θij) > 1e-15

            #output("    => Orbitals $i and $j to be rotated by $(θij)")
            # Rotate
            Cnew[:, i] =  cos(θij)*Cold[:,i] + sin(θij)*Cold[:,j]
            Cnew[:, j] = -sin(θij)*Cold[:,i] + cos(θij)*Cold[:,j]

            θvals[:,j] .= 0.0
            θvals[:,i] .= 0.0
            θvals[j,:] .= 0.0
            θvals[i,:] .= 0.0

            nrot += 1

            _, cmax = findmax(abs, θvals)
            θij = θvals[cmax]
            i,j = cmax.I
        end

        metric = pm_metric(Cnew, ints)
        δC = √(sum((Cnew .- Cold) .^2)/length(Cold))
        output("{:5d}     {:<10.5f}     {:<10.5f}     {:>5d}", iter, metric, δC, nrot)

        Cold .= Cnew
        iter += 1
    end
end

function get_Pst(C, ints)

    Nocc = size(C, 2)

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

    S = ints["S"]
    Qija = zeros(Nocc, Nocc, natoms)
    off = 0
    for a in 1:natoms
        n = basis_per_atom[a]
        for i in 1:Nocc
            for j in i:Nocc
                Ci = C[:,i]
                Cj = C[:,j]
                q = diag(Ci*transpose(Cj)*S + Cj*transpose(Ci)*S)
                Qija[i,j,a] = 0.5*sum(q[(off+1):(off+n)])
                if i != j
                    Qija[j,i,a] = Qija[i,j,a]
                end
            end
        end
        off += n
    end

    return Qija
end

function jacobi_sweep2(C, ints)

    Norbs = size(C, 2)
    Natoms = length(ints.molecule.atoms)

    output("Initiating Pipek-Mezey Localization")
    output(" • Algorithm: Jacobi Sweeps (Pair-wise rotations)\n")

    output("Initial PM metric: $(pm_metric(C, ints))")

    θvals = zeros(Norbs, Norbs)

    max_iter = 100
    iter = 1
    δC = 1.0
    while δC > 1e-10

        if iter > max_iter
            break
        end
        output("\n\n Iter $iter")

        # Get ⟨s|Pₐ|t⟩
        P = get_Pst(C, ints)

        for s = 1:Norbs
            for t = 1:(s-1)
                Ast = 0.0
                Bst = 0.0
                for a = 1:Natoms
                    Ast += P[s,t,a]^2 - 0.25*(P[s,s,a] - P[t,t,a])^2
                    Bst += P[s,t,a]*(P[s,s,a] - P[t,t,a])
                end
                AB = Ast^2 + Bst^2
                if abs(AB) > 1e-16
                    cos4θ = -Ast / √AB
                    sin4θ = Bst / √AB
                    #θvals[s,t] = 0.25 * acos(cos4θ) #* (Bst > 0 ? 1 : -1)
                    θvals[s,t] = 0.25*asin(sin4θ) + π/2
                else
                    θvals[s,t] = 0.0
                end
                #if s != t 
                #    θvals[t,s] = 0.0
                #end
            end
        end

        θmax, cmax = findmax(θvals)
        θmin, cmin = findmin(θvals)
        θ = abs(θmax) > abs(θmin) ? θmax : θmin
        i,j = abs(θmax) > abs(θmin) ? cmax.I : cmin.I

        output("Orbitals $i and $j to be rotated by $(θ)")

        # Rotate
        Cnew = deepcopy(C)
        Cnew[:, i] =  cos(θ)*C[:,i] + sin(θ)*C[:,j]
        Cnew[:, j] = -sin(θ)*C[:,i] + cos(θ)*C[:,j]

        output("Final PM metric: $(pm_metric(Cnew, ints))")
        δC = √(sum((Cnew .- C) .^2)/length(C))
        println(δC)

        C .= Cnew

        iter += 1
    end
    θvals
end