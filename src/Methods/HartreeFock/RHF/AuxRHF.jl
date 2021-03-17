# Auxiliar functions for RHF computations

"""
    RHFEnergy(D::FermiMDArray{Float64}, H::FermiMDArray{Float64}, F::FermiMDArray{Float64})

Given a density matrix `D`, Core Hamiltonian `H` and Fock matrix `F`, computes the RHF energy
"""
function RHFEnergy(D::FermiMDArray{Float64}, H::FermiMDArray{Float64},F::FermiMDArray{Float64})
    return sum(D .* (H .+ F))
end

"""
    build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ERI::FermiMDArray{Float64,4})

Builds a Fock matrix (into `F`) using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`.
"""
function build_fock!(F::FermiMDArray{Float64}, H::FermiMDArray{Float64}, D::FermiMDArray{Float64}, ERI::FermiMDArray{Float64,4})
    F .= H
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
end

"""
    build_fock!(F::FermiMDArray{Float64,2}, H::FermiMDArray{Float64,2}, D::FermiMDArray{Float64,2}, b::FermiMDArray{Float64,3})

Builds a Fock matrix (into `F`) using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`
approximated by density fitting.
"""
function build_fock!(F::FermiMDArray{Float64,2}, H::FermiMDArray{Float64,2}, D::FermiMDArray{Float64,2}, b::FermiMDArray{Float64,3})
    F .= H
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
end