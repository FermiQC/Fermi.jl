abstract type AbstractAtom end
begin
    struct Atom <: AbstractAtom
        Z::Int8
        position::Array{Float64,1}
    end
end
