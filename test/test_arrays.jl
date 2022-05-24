using LinearAlgebra
using TensorOperations
import Strided: UnsafeStridedView

@testset "Arrays" begin

    idx = [ (1,1,1,1), (2,2,2,2), (3,3,3,3)]
    val = rand(3)
    A = FermiSparse(idx, val)
    @test ndims(A) == 4
    @test length(A) == 3
    @test eltype(A) == Float64

end