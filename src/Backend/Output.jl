using Fermi.Options
using Formatting

export output

"""
    Fermi.output(str, x...; ending="\\n")

Writes a Python-style String `str` formatted with the values `x...`. The formatting is done through the package Formatting.jl.
The string is printed into the standard output or written into a file depending on the `printstyle` option. See below:

    printstyle: 
    "repl"  Default - Returns the string into the Julia REPL (standard output)
    "file"    Write the string into a file. Path must be specified with the keyword `output`
    "both"    Returns the value into stdout and also writes to file
    "none"    Does not write or print anything

In order to write to a file, the keyword `output` must be the path to the destination file.
The default path is "output.jl". For help setting options see `@set`.

# Examples
```
julia> @set {
    printstyle file
    output hydrogen.jl
}
julia> output("The Hydrogen Atom")
julia> output("1s Energy", ending=":") 
julia> output(" {:2.2f} hartrees",-0.5)
shell> cat hydrogen.jl
The Hydrogen Atom
1s Energy: -0.50 hartrees
```
"""
function output(str, x...; ending="\n")

    style = Options.get("printstyle")
    f = format(str*ending, x...)

    if style == "none"
        nothing
    elseif style == "repl"
        print(f)
    elseif style == "file"
        path = Options.get("output")
        open(path, "a") do io
            write(io, f)
            flush(io)
        end
    elseif style == "both"
        path = Options.get("output")
        open(path, "a") do io
            write(io, f)
            flush(io)
        end
        print(f)
    else
        throw(FermiException("printing style not recognized: $style. Accepted `printstyle` values: repl, file, both, and none"))
    end
end
