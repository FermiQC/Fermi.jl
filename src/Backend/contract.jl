function exclude(a::String,b::String)
    a2 = deepcopy(a)
    for idx in 1:length(b)
        i = b[idx:idx]
        if occursin(i,a2)
            a2 = replace(a2,i=>"")
        end
    end
    a2
end
function getdims(iC,iA,iB,A,B)
    sza = size(A)
    szb = size(B)
    szc = []
    for i in 1:length(iA)
        if occursin(iA[i:i],iC)
            push!(szc,sza[i])
        end
    end
    for i in 1:length(iB)
        if occursin(iB[i:i],iC)
            push!(szc,szb[i])
        end
    end
    szc
end
function contract(A,B,iA,iB)
    # figure out size of C, iC
    contract(A,B,1.0,1.0,iA,iB)
end
function contract(A,B,sA,sB,iA,iB)
    iC = exclude(iA,iB)*exclude(iB,iA)
    C = zeros(eltype(A),getdims(iC,iA,iB,A,B)...)
    contract!(C,A,B,1.0,sA,sB,iC,iA,iB)
    C
end


function contract!(C,A,B,iC,iA,iB)
    contract!(C,A,B,1.0,1.0,1.0,iC,iA,iB)
end

function contract!(C,A,B,sC::T,sA::T,sB::T,iC,iA,iB) where T <: AbstractFloat
    contract!(C,A,B,sC,sA,sB,iC,iA,iB,
              Fermi.ComputeEnvironment.interconnect,Fermi.ComputeEnvironment.communicator,
              Fermi.ComputeEnvironment.accelerator,Fermi.ComputeEnvironment.contractor)
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

function contract!(C::T1,
                   A::T2,
                   B::T3,  #tensors
                   sC,sA,sB, #scaling factors
                   iC,iA,iB, #indices
                   ic::I where I <: Fermi.Environments.AbstractInterconnect,
                   comm::Fermi.Environments.NoCommunicator,
                   acc::Fermi.Environments.NoAccelerator,
                   contractor::Fermi.Environments.TBLIS) where { T1 <: Union{Array,Fermi.MemTensor},
                                                                 T2 <: Union{Array,Fermi.MemTensor},
                                                                 T3 <: Union{Array,Fermi.MemTensor}}

    Ap = Fermi.data(A)
    Bp = Fermi.data(B)
    Cp = Fermi.data(C)
    _A = TBLIS.TTensor{eltype(Ap)}(Ap,sA)
    _B = TBLIS.TTensor{eltype(Bp)}(Bp,sB)
    _C = TBLIS.TTensor{eltype(Cp)}(Cp,sC)
    TBLIS.mul!(_C,_A,_B,iA,iB,iC)
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
