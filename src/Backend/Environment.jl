module Environments
import MPI

abstract type  AbstractInterconnect end
begin
    struct SSH_IC  <: AbstractInterconnect end
    struct MPI_IC  <: AbstractInterconnect end
    struct TCP_IC  <: AbstractInterconnect end
    struct No_IC   <: AbstractInterconnect end
end

abstract type AbstractCommunicator end
begin
    struct NoCommunicator <: AbstractCommunicator end
    struct MPI_World <: AbstractCommunicator
        comm::MPI.Comm
    end
    struct MachineList <: AbstractCommunicator
        machines::Array{String}
    end
    struct TCPList <: AbstractCommunicator
        ports::Array{String}
    end
end

abstract type AbstractAccelerator end
begin
    struct NoAccelerator <: AbstractAccelerator end
    abstract type AbstractGPU <: AbstractAccelerator end
    begin
        struct AMD_GPU <: AbstractGPU 
            memory::UInt128
        end
        struct NVIDIA_GPU <: AbstractGPU 
            memory::UInt128
        end
    end
    abstract type AbstractFPGA <: AbstractAccelerator end
    abstract type AbstractCoprocessor <: AbstractAccelerator end
end

abstract type AbstractContractor end
begin
    struct TTGT <: AbstractContractor end
    struct TBLIS <: AbstractContractor end
end

struct Environment
    interconnect::T where T <: AbstractInterconnect
    homogenous::Bool
    accelerator::T  where T <: AbstractAccelerator
    communicator::T         where T <: AbstractCommunicator
    contractor::T   where T <: AbstractContractor
end

function detect_Environment()
end
end #module
