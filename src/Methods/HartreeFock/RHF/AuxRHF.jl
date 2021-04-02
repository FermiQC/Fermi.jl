# Auxiliar functions for RHF computations
"""
    function RHF_core_guess(ints::IntegralHelper)

Compute an initial orbital coefficient matrix using the core guess.
"""
function RHF_core_guess(ints::IntegralHelper)
    output("Using Core Guess")
    S = ints["S"]
    Λ = S^(-1/2)
    F = ints["T"] + ints["V"]
    Ft = Λ*F*Λ'

    # Get orbital energies and transformed coefficients
    _, Ct = diagonalize(Ft, hermitian=true)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ
end

"""
    function RHF_gwh_guess(ints::IntegralHelper)

Compute an initial orbital coefficient matrix using the GWH guess.
"""
function RHF_gwh_guess(ints::IntegralHelper)

    # Form GWH guess
    output("Using GWH Guess")
    molecule = ints.molecule
    S = ints["S"]
    d, U = diagonalize(S, sortby = x->1/abs(x))
    Λ = S^(-1/2)
    idxs = [abs(d[i]) > 1E-7 for i = eachindex(d)]

    H = real.(ints["T"] + ints["V"])
    ndocc = molecule.Nα
    nvir = size(S,1) - ndocc
    F = similar(S)

    for i = 1:ndocc+nvir
        F[i,i] = H[i,i]
        for j = i+1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ'*F*Λ

    # Get orbital energies and transformed coefficients
    _, Ct = diagonalize(Ft, hermitian=true)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ         
end

"""
    RHFEnergy(D::FermiMDArray{Float64}, H::FermiMDArray{Float64}, F::FermiMDArray{Float64})

Compute RHF energy given a density matrix `D`, Core Hamiltonian `H` and Fock matrix `F`.
"""
function RHFEnergy(D::FermiMDArray{Float64}, H::FermiMDArray{Float64}, F::FermiMDArray{Float64})
    return sum(D .* (H .+ F))
end

"""
    build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ERI::FermiMDArray{Float64,4})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`.
"""
function build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ERI::FermiMDArray{Float64,4})
    F .= H
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
end

"""
    build_fock!(F::FermiMDArray{Float64,2}, H::FermiMDArray{Float64,2}, D::FermiMDArray{Float64,2}, b::FermiMDArray{Float64,3})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`
approximated by density fitting.
"""
function build_fock!(F::FermiMDArray{Float64,2}, H::FermiMDArray{Float64,2}, D::FermiMDArray{Float64,2}, b::FermiMDArray{Float64,3})
    F .= H
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
end
