abstract type AbstractReferenceWavefunction end
begin 
    mutable struct ReferenceWavefunction{T} <: AbstractReferenceWavefunction where T <: AbstractFloat
        refEnergy::T
        vnuc::T
        nocca::Int
        noccb::Int
        nvira::Int
        nvirb::Int
        basis::B where B <: Lints.BasisSet
        molecule::M where M <: Lints.Molecule
        Ca::Array{T,2} #AO->MO coefficients a
        Cb::Array{T,2} #AO->MO coefficients b
        S::Array{T,2}
        T::Array{T,2}
        V::Array{T,2}
        epsa::Array{T,1} #orbital eigenvalues a
        epsb::Array{T,1} #orbital eigenvalues b
        ERI::I where I <: Fermi.AbstractTensor#AO basis electron repulsion integrals
    end
end

function ReferenceWavefunction(basis,molecule,nocca,noccb)
    S_engine = Lints.OverlapEngine(nprim,l)
    T_engine = Lints.KineticEngine(nprim,l)
    V_engine = Lints.NuclearEngine(nprim,l,molecule)
    S = zeros(sz,sz)
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    Lints.make_2D(S,S_engine,basis)
    Lints.make_2D(T,T_engine,basis)
    Lints.make_2D(V,V_engine,basis)


