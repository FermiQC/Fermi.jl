module Error

export InvalidFermiOption
export MethodArgument
export InvalidMolecule

struct InvalidFermiOption <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidFermiOption) = print(io, "InvalidFermiOption: ", e.msg)

struct MethodArgument <: Exception
    msg::String
end
Base.showerror(io::IO, e::MethodArgument) = print(io, "MethodArgument: ", e.msg)

struct InvalidMolecule <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidMolecule) = print(io, "InvalidMolecule: ", e.msg)

end # module
