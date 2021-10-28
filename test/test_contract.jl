using Fermi.Options
using TensorOperations
import TBLIS
import TensorOperations: contract!

@testset "Tensor Contraction" begin
    Fermi.tblis_set_num_threads(1)
    @test TBLIS.get_num_threads() == 1

    # Using Aᵢₘₙ⋅Bₘₙⱼ = Cᵢⱼ
    oindA = (1,)
    cindA = (2,3)
    oindB = (3,)
    cindB = (1,2)
    tindC = (1,2)
    x = Fermi.oind2eins(oindA, cindA, oindB, cindB, tindC)
    @test x == ("Bbc", "bcC", "BC")

    A = FermiMDrand(2,5,5)
    B = FermiMDrand(5,5,2)
    C1 = FermiMDzeros(2,2)
    C2 = FermiMDzeros(2,2)

    @set tblis true
    α = rand()
    β = rand()
    contract!(α, A, :N, B, :N, β, C1, oindA, cindA, oindB, cindB, tindC)
    contract!(α, A.data, :N, B.data, :N, β, C2.data, oindA, cindA, oindB, cindB, tindC)
    @test C1 ≈ C2

    @test_throws FermiException contract!(α, A, :C, B, :N, β, C1, oindA, cindA, oindB, cindB, tindC)

    @set tblis false
    C2 = FermiMDzeros(2,2)
    contract!(α, A, :N, B, :N, β, C2, oindA, cindA, oindB, cindB, tindC)
    @test C1 ≈ C2

    A = FermiMDrand(5,5)
    B = FermiMDrand(5,5)
    C = zeros(5,5)

    @set tblis true
    @tensor C[i,j] = A[i,k]*B[k,j]
    @test C ≈ A*B

    A = rand(5,5)
    B = FermiMDrand(5,5)
    C = zeros(5,5)

    @tensor C[i,j] = A[i,k]*B[k,j]
    @test C ≈ A*B

    @tensor C[i,j] = B[i,k]*A[k,j]
    @test C ≈ B*A
end
