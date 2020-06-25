module Input

using Fermi
using Fermi.HartreeFock.RHF
using Lints

export run #shorthand for read->exec
export read
export exec

include("psi4keywords.dat")

"""
    run(fname::String)

read and execute the input file specified by fname. 
"""
function run(fname::String = "input.dat")
    cont = read(fname)
    exec(cont)
end

"""
    read(fname::String)

read an input file and convert to a Dict. exec takes the output of this function as input.
"""
function read(fname)
    if isfile(fname)
        lines = readlines(fname)
    else
        error("Input file not found $fname")
    end
    molecule,options,command = parse(lines)
    cont = Dict("molecule"=>molecule,"options"=>options,"command"=>command)
end

"""
    exec(cont::Dict)

execute the computation specified in cont.
cont is composed of:
    Dict("molecule"=>molecule::String,
         "options"=>options::Dict,
         "command"=>command::function)
"""
function exec(cont)
    open("/tmp/molfile.xyz","w") do molfile
        write(molfile,"3\n\n")
        write(molfile,cont["molecule"])
    end
    Lints.libint2_init()
    mol = Lints.Molecule("/tmp/molfile.xyz")
    bas = Lints.BasisSet(cont["options"]["basis"],mol)
    Fermi_options = Dict{Any,Any}(
                       :doprint=>false,
                       :maxit=>40,
                       :return_T=>false,
                       :quiet=>true,
                       :frozen => 1,
                       :active => 8
                      )
    rhfwfn = Fermi.ReferenceWavefunction(bas,mol,5,5)
    do_RHF(rhfwfn)
    com = cont["command"]
    E = com(rhfwfn; Fermi_options...)
    #Lints.libint2_finalize()
    return E
end

function parse(lines;delim="--")
    chunks = split(join(lines,"\n"),delim)
    for chunk in chunks
        expr = Meta.parse(chunk)
        eval(expr)
    end
    molecule,options,command
end


end # module
