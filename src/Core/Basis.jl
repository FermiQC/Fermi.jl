abstract type AbstractBasis end
struct Basis <: AbstractBasis
    lbasis::B where B <: Lints.BasisSet
    bname::String
end
struct NullBasis <: AbstractBasis end
