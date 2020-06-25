function contract(A,B,iA,iB)
    # figure out size of C, iC
    # C = 
    # iC = 
    contract!(C,A,B,iC,iA,iB)
end

function contract!(C,A,B,iC,iA,iB)
    contract!(C,A,B,1.0,1.0,1.0,iC,iA,iB)
end

function contract!(C,A,B,sC::T,sA::T,sB::T,iC,iA,iB) where T <: AbstractFloat
    contract!(C,A,B,sC,sA,sB,iC,iA,iB,Fermi.Env.interconnect,Fermi.Env.communicator,Fermi.Env.contractor)
    return C
end

function contract!(C,A,B, #tensors
                   sC,sA,sB, #scaling factors
                   iC,iA,iB, #indices
                   ic::I where I <: Fermi.Environments.AbstractInterconnect,
                   comm::Fermi.Environments.NoCommunicator,
                   acc::Fermi.Environments.NoAccelerator,
                   contractor::Fermi.Environments.TTGT)
    # transpose
    # transpose
    # GEMM
    # transpose
end

function contract!(C,A,B, #tensors
                   sC,sA,sB, #scaling factors
                   iC,iA,iB, #indices
                   ic::I where I <: Fermi.Environments.AbstractInterconnect,
                   comm::Fermi.Environments.NoCommunicator,
                   acc::Fermi.Environments.NoAccelerator,
                   contractor::Fermi.Environments.TBLIS)
    mul!(C,A,B,iA,iB,iC)
end

function contract!(C,A,B, #tensors
                   sC,sA,sB, #scaling factors
                   iC,iA,iB, #indices
                   ic::Fermi.Environments.MPI_IC,
                   comm::Fermi.Environments.MPI_World,
                   acc::Fermi.Environments.AMD_GPU,
                   contractor::Fermi.Environments.TTGT)
    # Super special MPI aware GPU code specialized for AMD products
end
