module Integrals

using Fermi
using Lints
using LinearAlgebra

abstract type AbstractIntegrals end

struct AOIntegrals{T} <: AbstractIntegrals where T <: AbstractFloat
    Vnuc::T
    S::Array{T,2}
    T::Array{T,2}
    V::Array{T,2}
    ERI:: I where I <: Fermi.AbstractTensor
end


function AOIntegrals()

    mol = Fermi.CurrentOptions["molstring"]
    basis = Fermi.CurrentOptions["basis"]
    AOIntegrals(mol, basis)
end

function AOIntegrals(mol::String, basis::String)

    charge = Fermi.CurrentOptions["charge"]
    multiplicity = Fermi.CurrentOptions["multiplicity"]
    noccα, noccβ = get_num_electrons(mol, charge, multiplicity)
    AOIntegrals(mol, basis, noccα, noccβ)
end

function AOIntegrals(mol::String, basis::String, noccα::Int, noccβ::Int)

    AOIntegrals(basis, mol, noccα, noccβ,
                          Fermi.ComputeEnvironment.interconnect,
                          Fermi.ComputeEnvironment.communicator,
                          Fermi.ComputeEnvironment.accelerator)
end

function AOIntegrals(mol::String, basis::String, noccα::Int, noccβ::Int, unit::String,
                               interconnect::Fermi.Environments.No_IC,
                               communicator::Fermi.Environments.NoCommunicator,
                               accelerator::Fermi.Environments.NoAccelerator)

    Vnuc = nuclear_repulsion(mol, unit)

    if unit == "bohr"
        molecule = convert_xyz_unit(mol, from="bohr", to="angstrom")
    else
        molecule = mol
    end

    open("/tmp/molfile.xyz","w") do molfile
        natom = length(split(strip(molecule),"\n"))
        write(molfile,"$natom\n\n")
        write(molfile,molecule)
    end

    Lints.libint2_init()
    mol = Lints.Molecule("/tmp/molfile.xyz")
    bas = Lints.BasisSet(basis, mol)

    nprim = Lints.max_nprim(bas)
    l = Lints.max_l(bas)
    S_engine = Lints.OverlapEngine(nprim,l)
    T_engine = Lints.KineticEngine(nprim,l)
    V_engine = Lints.NuclearEngine(nprim,l,mol)
    I_engines = []
    sz = Lints.getsize(S_engine,bas)
    for i in 1:Threads.nthreads()
        push!(I_engines,Lints.ERIEngine(nprim,l))
    end
    S = zeros(sz,sz)
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    I = zeros(sz,sz,sz,sz)
    Lints.make_2D(S,S_engine,bas)
    Lints.make_2D(T,T_engine,bas)
    Lints.make_2D(V,V_engine,bas)
    Lints.make_ERI(I,I_engines,bas)
    I = Fermi.MemTensor(I)
    Lints.libint2_finalize()

    return AOIntegrals{Float64}(Vnuc, S, T, V, I)
end

end #module
