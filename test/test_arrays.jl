using LinearAlgebra
using TensorOperations
import Strided: UnsafeStridedView

@testset "Arrays" begin

    A = FermiMDrand(5,5)
    @test FermiMDArray(A) === A
    @test size(A) == (5,5)
    @test A[2,2] == A.data[2,2]
    A[2,2] = 0.75
    @test A[2,2] == 0.75
    @test ndims(A) == 2
    @test length(A) == 25
    At = similar(A)
    At .= transpose(A)
    @test permutedims(A, (2,1)) == At
    permutedims!(A, (2,1)) 
    @test A == At

    @test diag(A) == diag(A.data)

    # Make it Hermitian
    A = A + transpose(A)

    @test A^2 == A.data^2
    @test A^-1 ≈ A^-1.0

    Hvals, Hvec = diagonalize(Hermitian(A), hermitian=true)
    adjoint(Hvec) * A*Hvec ≈ diagm(Hvals)
    vals, vec = diagonalize(A, hermitian=false)
    vec^-1 * A*vec ≈ diagm(vals)

    b = rand(5,5)
    B = FermiMDArray(b)
    @test A + b == A + B == b + A
    @test A - b == A - B == -(b - A)
    @test 2*A == 2*A.data == A*2

    @test FermiMDArray(2.0) == 2.0

    C = FermiMDzeros(5,5)
    @test C == zeros(5,5)

    Bv = similar(B, (5,))
    Bv .= B[1,:]
    @test Bv == B.data[1,:]

    Ct = rand(Float32, (5,5))
    C .= Ct
    @test C.data == Base.convert(Array{Float64}, Ct)
    
    @test UnsafeStridedView(A) == UnsafeStridedView(A.data)

    idx = [ (1,1,1,1), (2,2,2,2), (3,3,3,3)]
    val = rand(3)
    A = FermiSparse(idx, val)
    @test ndims(A) == 4

    # x = @capture_out display(A)
    # println("")
    # @test occursin("Fermi Memory-held Dense Array", x)
end