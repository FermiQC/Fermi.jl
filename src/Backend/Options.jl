export @get
export @set
export @reset
export @molecule
export @lookup

using Formatting
using PrettyTables

"""
    Fermi.Options

Module to manage options in Fermi. 

# Functions

    Fermi.set(option, value)     Set an <option> to a given <value>
    Fermi.get(option)            Return the current value of an <option>
    Fermi.reset()                Reset all options to default values
    Fermi.reset(option)          Reset a specific <option> to its default value
    Fermi.molecule(molstring)    Read in a String for the `molstring` option

Alternatively, at global scope, one can use the corresponding macros that create shortcuts
for the commands above

# Macros

    @set <option> <value>   Set an <option> to a given <value>
    @get <option>           Return the current value of an <option>
    @reset                  Reset all options to default values
    @reset <option>         Reset a specific <option> to its default value
    @molecule               Read in a String for the `molstring` option
"""
module Options
using Fermi.Error

"""
    Fermi.Options.Default

Dictionary containing default options for Fermi. Any information not given
explicitly to Methods is obtained from here.
"""
const Default = Dict{String,Union{Float64,Int,String,Bool}}(
                                  "molstring" => """
                                  O        1.2091536548      1.7664118189     -0.0171613972
                                  H        2.1984800075      1.7977100627      0.0121161719
                                  H        0.9197881882      2.4580185570      0.6297938832
                                  """,
                                  "printstyle" => "repl",
                                  "output" => "output.jl",
                                  "basis" => "sto-3g",
                                  "jkfit" => "auto",
                                  "rifit" => "auto",
                                  "charge" => 0,
                                  "multiplicity" => 1,
                                  "unit" => "angstrom",
                                  "reference" => "rhf",
                                  "scf_max_iter" => 50,
                                  "scf_max_rms" => 10^-9,
                                  "scf_e_conv" => 10^-10,
                                  "scf_alg" => 1,
                                  "mp2_alg" => 1,
                                  "df" => false, 
                                  "oda" => true,
                                  "oda_cutoff" => 1E-1,
                                  "oda_shutoff" => 20,
                                  "scf_guess" => "gwh",
                                  "quiet" => true,
                                  "precision" => "double",
                                  "cc_alg" => 1,
                                  "ci_alg" => "aci",
                                  "det_size" => 64,
                                  "σ"      => 0.001,
                                  "γ"      => 1.0,
                                  "ζ"      => 0.95,
                                  "ζsize"  => 0,
                                  "d_conv" => 8,
                                  "cc_max_iter" => 50,
                                  "cc_max_rms" => 10^-10,
                                  "cc_e_conv" => 10^-10,
                                  "bcc_max_t1" => 1^-7,
                                  "drop_occ" => 0,
                                  "drop_vir" => 0,
                                  "diis" => true,
                                  "cc_diis" => true,
                                  "ndiis" => 8,
                                  "cc_ndiis" => 3,
                                  "bcc_tol" => 1.0e-6,
                                  "diis_prec" => "single",
                                  "diis_start" => 3,
                                  "cc_damp_ratio" => 0.0,
                                  "pt_alg" => 1,
                                  "num_frozen" => 0,
                                  "aci_print_screen" => 0,
                                  "cas_frozen" => 0,
                                  "cas_active" => -1,
                                  "cas_cutoff" => 10^-9,
                                  "cas_nroot" => 1,
                                  "min_matrix_elem" => 10^-9,
                                  "precision_override" => false,
                                  "tblis" => false
                                 )
"""
    Fermi.Options.Current

Dictionary containing user defined options for Fermi. Unspecified options are obtained from 
Fermi.Default.
"""
Current = Dict{String,Union{Float64,Int,String,Bool}}()


"""
    Fermi.Options.get(key)

Return the current options set for the given `key`.
See also `@get`
"""
function get(key::String)

    key = lowercase(key)
    if haskey(Current, key)
        return Current[key]
    elseif haskey(Default, key)
        return Default[key]
    else
        throw(InvalidFermiOption(key*" is not a valid option."))
    end
end


"""
    Fermi.Options.set(key, val)

Set the option `key` to the value `val`.
See also `@set`
"""
function set(key::String, val::Union{String, Bool, Float64, Int})
    # Make all string lowercase
    key = lowercase(key)
    if key != "molstring" && val isa String
        val = lowercase(val)
    end
    if !(haskey(Default, key))
        throw(InvalidFermiOption(key*" is not a valid option."))
    end

    # Check if the type of the variable passed matches the type
    # saved in the Default dictionary
    dtype = typeof(Default[key])
    if !isa(val, dtype)
        # If it does not match...

        # If the expected type is not a String we try to make a conversion
        if dtype !== String
            try newval = dtype(val)
                # Recall the function with the new variable with adjusted type
                return set(key, newval)
            catch InexactError
                # If the conversion is not possible
                throw(InvalidFermiOption("Impossible to convert data type for $key. Expected: $dtype, Got: $(typeof(val))"))
            end
        end
        # Otherwise, throw an error
        throw(InvalidFermiOption("Invalid data type for $key. Expected: $dtype, Got: $(typeof(val))"))
    else
        Current[key] = val
    end
    get(key)
end

"""
    Fermi.Options.reset(key...)

Reset all given `key` to the default values.
See also `@reset`
"""
function reset(key::String="all")
    if key == "all"
        for k in keys(Current)
            delete!(Current, k)
        end
    elseif haskey(Current,key)
        delete!(Current, key)
    elseif !haskey(Default, key)
        throw(InvalidFermiOption(key*" is not a valid option."))
    end
end
end #module

# Aux function to convert Symbols to Strings
function symbol_to_string(S)
    S = repr(S)
    return String(strip(filter(c->!occursin(c," {}():"),S)))
end

"""
    Fermi.@get

Returns the current value of an option. 

# Examples

```
julia> @get basis
"sto-3g"      # Default

julia> @set basis cc-pVDZ
julia> @get basis
"cc-pVDZ"
```
Without arguments, it returns all user defined keywords
```    
julia> @set {
    basis cc-pVDZ
    scf_max_rms 10^-8
    diis true
}
julia> @get
┌─────────────┬───────────────┐
│     Keyword │ Current Value │
├─────────────┼───────────────┤
│ scf_max_rms │   1.00000e-08 │
│        diis │          true │
│       basis │       cc-pvdz │
└─────────────┴───────────────┘
```
"""
macro get(opt)
    A = symbol_to_string(opt)
    quote
        Fermi.Options.get($A)
    end |> esc
end

macro get()
    quote
        Keys = collect(keys(Fermi.Options.Current))
        l = length(Keys)
        if l == 0
            println("No user defined keywords found.")
        else
            data = Array{Union{Nothing,String},2}(nothing,l,2)
            for i in eachindex(Keys)
                k = Keys[i]
                data[i,1] = k
                val = Fermi.Options.get(k) 
                if val isa Bool || val isa Int
                    data[i,2] = "$val"
                elseif val isa String
                    # If multiple lines (e.g. molstring) display only the first line
                    if occursin("\n", val)
                        lines = split(val, "\n")
                        data[i,2] = lines[1] * " [+$(length(lines)-1) lines]"
                    else 
                        data[i,2] = "$val"
                    end
                else
                    data[i,2] = format("{:5.5e}", val)
                end
            end
            pretty_table(data; header=["Keyword", "Current Value"])
        end
    end
end

"""
    Fermi.@set

Set options for Fermi computations. It saves the options into Fermi.CurrentOptions.

*Usage:*  @set A B

A is set to B. By default A is taken as a string. 

# Examples

```
@set basis cc-pVDZ
@set cc_max_iter 100
@set e_conv 10^-9
```
One can also use the block syntax
```
@set {
    basis cc-pVDZ
    cc_max_iter 100
    e_conv 10^-9
}
```
Note for basis set: Having * at the end of a line is considered an incomplete expression by Julia. Thus, for basis such as 6-31g*
you should use quotes
```
@set basis "6-31g"
```
"""
macro set(A::Symbol, B::Symbol)
    clean_up(s) = String(filter(c->!occursin(c," ():"),s))
    key = clean_up(repr(A))
    val = clean_up(repr(B))
    return quote
        Fermi.Options.set($key,$val)
    end |> esc
end

macro set(A::Symbol, B::Union{Number, Bool, String})
    clean_up(s) = String(filter(c->!occursin(c," ():"),s))
    key = clean_up(repr(A))
    quote
        Fermi.Options.set($key, $B)
    end
end

macro set(A::Symbol, B::Expr)
    clean_up(s) = String(filter(c->!occursin(c," ():"),s))
    key = clean_up(repr(A))

    if all(x->typeof(x)<:Number, B.args[2:end])
        quote
            Fermi.Options.set($key, $B)
        end
    else
        val = clean_up(repr(B))
        quote
            Fermi.Options.set($key, $val)
        end
    end
end

macro set(block::Expr)

    # Create an empty quote for the output
    out = quote end

    # Each argument in the expression block is a line
    # e.g. 
    # @set {
    #    A B   -> First line (Expr object)
    #    C D   -> Second line (Expr object)
    #}
    lines = block.args

    for line in lines
        # args represent the arguments of each line. In the example above
        # args = {:A, :B} and {:C, :D}
        args = line.args

        # We must assert there are only two arguments
        length(args) == 2 || throw(InvalidFermiOption("Too many arguments for @set: $line"))

        # A is the dictionary key, which is always taken as a String
        key = symbol_to_string(args[1])

        # B is the val
        val = args[2]
        valtype = typeof(val)

        # First we check the easy case: B is a Number or Bool
        if valtype <: Union{Number, Bool, String}
            push!(out.args, quote Fermi.Options.set($key, $val) end)

        # Now the case where B is a Symbol. 
        elseif valtype === Symbol
            # Turn it into a String
            val = symbol_to_string(val)
            push!(out.args, quote Fermi.Options.set($key,$val) end)

        # Finally, if the value is an Expression
        elseif valtype === Expr

            # If the expression is purely numerical, we just used it as it is
            if all(x->typeof(x)<:Number, val.args[2:end])
                push!(out.args, quote Fermi.Options.set($key,$val) end)

            # If not, String it out
            else
                val = symbol_to_string(val)
                push!(out.args, quote Fermi.Options.set($key,$val) end)
            end
        else
            # Outside those cases the user did something sketchy
            throw(InvalidFermiOption("Invalid input for $key: $val"))
        end
    end

    return out
end

"""
    Fermi.@reset(x="all", y...)

Restores any number of option to default setting. If one argument is "all", restores all options to their default values.

# Examples
```
julia> @set basis 6-31g
julia> @get basis
"6-31g"
julia> @reset basis 
julia> @get basis
"sto-3g"

julia> @set {
            scf_guess core
            diis false
        }
julia> @reset
julia> @get scf_core
"gwh"
julia> @get diis
true
```
"""
macro reset(x="all", y...)
    Keys = [String(x)]
    for i in y
        push!(Keys, String(i))
    end

    quote
        for k in $Keys
            Fermi.Options.reset(k)
        end
    end
end

"""
    Fermi.@molecule

Set the molecule to be used in computations.

# Example
```
@molecule {
C                 -0.00000000     0.00000000    -0.00000000
O                  0.00000000     0.00000000    -2.19732722
O                 -0.00000000    -0.00000000     2.19732722
}
```
"""
macro molecule(block)
    clean_up(s) = strip(filter(c->!occursin(c,"{}():"),s))
    mol = repr(block)
    mol = replace(mol, ";"=>"\n")
    mol = String(clean_up(mol))
    quote
        Fermi.Options.set("molstring", $mol)
    end
end

"""
    Fermi.@lookup

Look up valid keywords for Fermi.

# Example
```
julia> @lookup conv
┌────────────┬───────────────┐
│    Keyword │ Current Value │
├────────────┼───────────────┤
│  cc_e_conv │   1.00000e-10 │
│     d_conv │             8 │
│ preconv_t1 │          true │
│ scf_e_conv │   1.00000e-01 │
└────────────┴───────────────┘
```
"""
macro lookup(A::Symbol)
    clean_up(s) = strip(filter(c->!occursin(c,"{}():"),s))
    A = clean_up(repr(A))
    quote
        Keys = collect(keys(Fermi.Options.Default))
        filter!(s->occursin($A, s), Keys)
        l = length(Keys)
        if l == 0
            println("No keywords found containing: $($A)")
        else
            data = Array{Union{Nothing,String},2}(nothing,l,2)
            for i in eachindex(Keys)
                k = Keys[i]
                data[i,1] = k
                val = Fermi.Options.get(k) 
                if val isa Bool || val isa Int
                    data[i,2] = "$val"
                elseif val isa String
                    # If multiple lines (e.g. molstring) display only the first line
                    if occursin("\n", val)
                        lines = split(val, "\n")
                        data[i,2] = lines[1] * " [+$(length(lines)-1) lines]"
                    else
                        data[i,2] = "$val"
                    end
                else
                    data[i,2] = format("{:5.5e}", val)
                end
            end
            pretty_table(data; header=["Keyword", "Current Value"])
        end
    end
end
