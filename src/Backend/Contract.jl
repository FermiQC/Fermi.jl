import TensorOperations: contract!
import TensorOperations: IndexTuple
import TensorOperations: IndexError
import TBLIS

function oind2eins(oindA::NTuple{NAo}, cindA::NTuple{NAc},
          oindB::NTuple{NBo}, cindB::NTuple{NBc},
          tindC::NTuple{NCt}) where {NAo, NAc, NBo, NBc, NCt}

    # Check contraction conssitency.
    NAo + NBo == NCt || throw(IndexError("number of outer index not consistent."))
    NAc == NBc || throw(IndexError("number of contracted index not consistent."))

    cPadding = 'a' - 'A'
    einA = zeros(Int8, NAo+NAc)
    einB = zeros(Int8, NBo+NBc)

    # Outer indices.
    for i = 1:NAo
        einA[oindA[i]] = i
    end
    for i = 1:NBo
        einB[oindB[i]] = i + NAo
    end

    # Contracted indices.
    for i = 1:NAc
        einA[cindA[i]] = i + cPadding
        einB[cindB[i]] = i + cPadding
    end

    einA = string((einA .+'A')...)
    einB = string((einB .+'A')...)
    einC = string((tindC.+'A')...) # C has direct conversion relations.
    einA, einB, einC
end

function contract!(α, 
          A::FermiMDArray{T}, conjA::Symbol,
          B::FermiMDArray{T}, conjB::Symbol,
          β, 
          C::FermiMDArray{T}, 
          oindA::IndexTuple, cindA::IndexTuple, 
          oindB::IndexTuple, cindB::IndexTuple,
          tindC::IndexTuple, idk::IndexTuple, syms::Union{Nothing, NTuple{3,Symbol}} = nothing) where{T<:AbstractFloat}

    if Fermi.CurrentOptions["tblis"]
        # Check permutation consistency.
        # This check is copied from stridedarray.jl
        pA = (oindA...,cindA...)
        (length(pA) == ndims(A) && isperm(pA)) ||
            throw(IndexError("invalid permutation of length $(ndims(A)): $pA"))
        pB = (oindB...,cindB...)
        (length(pB) == ndims(B) && isperm(pB)) ||
            throw(IndexError("invalid permutation of length $(ndims(B)): $pB"))
        (length(tindC) == ndims(C) && isperm(tindC)) ||
            throw(IndexError("invalid permutation of length $(ndims(C)): $tindC"))

        # Complex conjugate is not supported by TBLIS yet.
        # Here we do a intermediate copying-conversion.
        if conjA == :C || conjA == :A
            A = conj(A)
        end
        if conjB == :C || conjB == :A
            B = conj(B)
        end
        einA, einB, einC = oind2eins(oindA, cindA, 
                                     oindB, cindB, 
                                     tindC)

        CType = eltype(C.data)
        _A = TBLIS.TTensor{CType}(A.data,CType(α))
        _B = TBLIS.TTensor{CType}(B.data)
        _C = TBLIS.TTensor{CType}(C.data,CType(β))

        TBLIS.mul!(_C, _A, _B, einA, einB, einC)
        return C.data
    else
        contract!(α, A.data, conjA, B.data, conjB, β, C.data, oindA, cindA, oindB, cindB, tindC)
    end
end

