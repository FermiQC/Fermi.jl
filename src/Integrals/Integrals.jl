module Integrals

using Fermi
using Lints

abstract type AbstractIntegrals end

struct AOIntegrals{T} <: AbstractIntegrals where T <: AbstractFloat
    S::Array{T,2}
    T::Array{T,2}
    V::Array{T,2}
    ERI:: I where I <: Fermi.AbstractTensor
end

#function AOIntegrals()
#    
#    molecule = Fermi.CurrentOptions["molstring"]
#    basis = Fermi.CurrentOptions["basis"]
#     
#
#end

function AOIntegrals(basis::String, molecule::String, noccα::Int, noccβ::Int)
    AOIntegrals(basis, molecule, noccα, noccβ,
                          Fermi.ComputeEnvironment.interconnect,
                          Fermi.ComputeEnvironment.communicator,
                          Fermi.ComputeEnvironment.accelerator)
end

function AOIntegrals(basis::String, molecule::String, noccα::Int,noccβ::Int,
                               interconnect::Fermi.Environments.No_IC,
                               communicator::Fermi.Environments.NoCommunicator,
                               accelerator::Fermi.Environments.NoAccelerator)


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
    Ca = zeros(size(S))
    Cb = zeros(size(S))
    T = zeros(sz,sz)
    V = zeros(sz,sz)
    I = zeros(sz,sz,sz,sz)
    Lints.make_2D(S,S_engine,bas)
    Lints.make_2D(T,T_engine,bas)
    Lints.make_2D(V,V_engine,bas)
    Lints.make_ERI(I,I_engines,bas)
    I = Fermi.MemTensor(I)
    Lints.libint2_finalize()

    return AOIntegrals{Float64}(S, T, V, I)
end





end #module
