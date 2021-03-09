module Error

export InvalidFermiOption
export MethodArgument
export InvalidMolecule

"""
    InvalidFermiOption

Error flag used when an invalid option is passed to Fermi.
"""
struct InvalidFermiOption <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidFermiOption) = print(io, "InvalidFermiOption: ", e.msg)

"""
    MethodArgument

Error flag used when an invalid arguments are passed into a Fermi function (e.g. RHF)
"""
struct MethodArgument <: Exception
    msg::String
end
Base.showerror(io::IO, e::MethodArgument) = print(io, "MethodArgument: ", e.msg)

"""
    InvalidMolecule

Error flag used when it is not possible to construct a molecule with the 
information given.
"""
struct InvalidMolecule <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidMolecule) = print(io, "InvalidMolecule: ", e.msg)

end # module
