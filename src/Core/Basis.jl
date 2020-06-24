abstract type AbstractBasis end
begin
    struct Basis <: AbstractBasis
        lbasis::B <: Lints.BasisSet
        bname::String
    end
    struct NullBasis <: AbstractBasis end
end
