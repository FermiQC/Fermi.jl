module HartreeFock
using Fermi
using Lints
#using PyCall
using TensorOperations
using LinearAlgebra
using Fermi.Output

include("RHF.jl")

function print_header()
    @output repeat("=",80)*"\n"
    @output "|    {:<74}|\n" "Hartree Fock"
    @output "|        {:<70}|\n" "Module by M.M. Davis"
    @output repeat("=",80)*"\n"
end
defaults = Dict{Any,Any}(
                        :e_convergence => 1E-10,
                        :d_convergence => 1E-10
                        )
end #module

