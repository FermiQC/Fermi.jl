export @reset
export @set
export @get
export @molecule

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
"""
macro get(opt)
    clean_up(s) = String(filter(c->!occursin(c," ():"),s))
    A = clean_up(repr(opt))
    quote
        Fermi.Options.get($A)
    end |> esc
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
B is evaluated at runtime. If the evaluation is not possible, B is converted to a string.
# Examples
```
mybasis = "6-31g"
julia> @set basis mybasis
"6-31g"
julia> @set basis yourbasis
"yourbasis"
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
macro set(A, B)
    clean_up(s) = String(filter(c->!occursin(c," ():"),s))
    key = clean_up(repr(A))
    val = clean_up(repr(B))
    
    quote 
        try
            Fermi.Options.set($key, $(esc(Meta.parse(val))))
        catch UndefVarError
            Fermi.Options.set($key,$val)
        end
    end
end

macro set(block)
    lines = split(repr(block),";")
    clean_up(s) = String(strip(filter(c->!occursin(c," {}():"),s)))
    out = quote end
    for l in lines
        key, val = split(strip(l), " ", limit=2)
        key = clean_up(key)
        val = clean_up(val)

        y = quote
            try
                Fermi.Options.set($key, $(esc(Meta.parse(val))))
            catch UndefVarError
                Fermi.Options.set($key,$val)
            end
        end

        for com in y.args
            push!(out.args, com)
        end
    end

    return out
end

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
const Default = Dict{String,Union{Float64,Int,String,Bool,Nothing}}(
                                  "molstring" => """
                                  O        1.2091536548      1.7664118189     -0.0171613972
                                  H        2.1984800075      1.7977100627      0.0121161719
                                  H        0.9197881882      2.4580185570      0.6297938832
                                  """,
                                  "basis" => "sto-3g",
                                  "jkfit" => "auto",
                                  "rifit" => "auto",
                                  "charge" => 0,
                                  "multiplicity" => 1,
                                  "unit" => "angstrom",
                                  "reference" => "rhf",
                                  "scf_max_iter" => 50,
                                  "scf_max_rms" => 10^-10,
                                  "scf_type" => "conventional",
                                  "oda" => true,
                                  "oda_cutoff" => 1E-1,
                                  "oda_shutoff" => 20,
                                  "scf_guess" => "gwh",
                                  "quiet" => true,
                                  "mp2_type" => "df",
                                  "precision" => "double",
                                  "cc_alg" => "CTF",
                                  "ci_alg" => "aci",
                                  "det_size" => 64,
                                  "σ"      => 0.001,
                                  "γ"      => 1.0,
                                  "ζ"      => 0.95,
                                  "ζsize"  => nothing,
                                  "e_conv" => 10,
                                  "d_conv" => 8,
                                  "cc_max_iter" => 50,
                                  "cc_max_rms" => 10^-10,
                                  "cc_e_conv" => 10^-10,
                                  "bcc_max_t1" => 1^-7,
                                  "preconv_t1" => true,
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
                                  "cc_diis_relax" => 3,
                                  "num_frozen" => 0,
                                  "aci_print_screen" => nothing,
                                  "cas_frozen" => 0,
                                  "cas_active" => -1,
                                  "cas_cutoff" => 10^-9,
                                  "cas_nroot" => 1,
                                  "min_matrix_elem" => 10^-9,
                                  "precision_override" => false,
                                  "tblis" => false
                                 )
"""
    Fermi.Current

Dictionary containing user options for Fermi. Unspecified options are obtained from 
Fermi.Default.
"""
Current = Dict{String,Union{Float64,Int,String,Bool,Nothing}}()

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

function set(key::String, val::Union{String, Bool, Float64, Int, Nothing})
    key = lowercase(key)
    if !(haskey(Default, key))
        throw(InvalidFermiOption(key*" is not a valid option."))
    else
        Current[key] = val
    end
    get(key)
end

"""
    Fermi.@reset(key="all")

Restores `key` option to default setting. If `key == "all"`, restores all options to their default values.
"""
macro reset(key="all") 
    if key == "all"
        return quote
            for key in keys(Fermi.CurrentOptions)
                Fermi.CurrentOptions[key] = Fermi.DefaultOptions[key]
            end
        end
    elseif String(key) in keys(Fermi.CurrentOptions)
        return quote
            Fermi.CurrentOptions[$(String(key))] = Fermi.DefaultOptions[$(String(key))]
        end
    else
        return quote
            throw(InvalidFermiOption("Invalid option "*$(String(key))))
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
    mol = clean_up(mol)
    CurrentOptions["molstring"] = String(mol)
end
end #module
