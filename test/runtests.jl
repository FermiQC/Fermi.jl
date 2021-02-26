using Fermi
using Test

tol = 1E-8
#include("test_Options.jl")
@testset "Fermi" begin
    include("test_RHF.jl")
    #include("test_CC.jl")
    #include("test_MP.jl")
    #include("test_Options.jl")
    #include("test_CI.jl")
end
