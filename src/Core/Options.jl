# default options go here
# if user specifies option they will be overwritten
export CurrentOptions
export InvalidFermiOption
export @reset
export @set
export @get
export @molecule

const DefaultOptions = Dict{String,Union{Float64,Int,String,Bool,Nothing}}(
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
                                  "scf_alg" => "conventional",
                                  "oda" => true,
                                  "oda_cutoff" => 1E-1,
                                  "oda_shutoff" => 20,
                                  "scf_guess" => "gwh",
                                  "quiet" => true,
                                  "mp2_type" => "df",
                                  "precision" => "double",
                                  "cc_alg" => "CTF",
                                  "ci_alg" => "aci",
                                  "σ"      => 0.001,
                                  "γ"      => 1.0,
                                  "e_conv" => 10,
                                  "d_conv" => 8,
                                  "cc_max_iter" => 50,
                                  "cc_max_rms" => 10^-10,
                                  "cc_e_conv" => 10^-10,
                                  "bcc_max_t1" => 1^-7,
                                  "preconv_T1" => true,
                                  "drop_occ" => 0,
                                  "drop_vir" => 0,
                                  "diis" => true,
                                  "cc_diis" => true,
                                  "ndiis" => 8,
                                  "cc_ndiis" => 3,
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
                                  "precision_override" => false
                                 )
"""
    Fermi.CurrentOptions

Dictionary containing options for Fermi. Any information not given
explicitly to Methods is obtained from here.
"""
CurrentOptions = deepcopy(DefaultOptions)

"""
    Fermi.@reset

Restores all options to their default values.
"""
macro reset(key=nothing) 
    if key == nothing
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
            error("Invalid option "*$(String(key)))
        end
    end
end



"""
    Fermi.@set

Set options for Fermi computations. It saves the options into Fermi.CurrentOptions.

*Usage:*  @set A B

A is set to B. By default A is taken as a string. B is evaluated at runtime. If the evaluation
is not possible, B is converted to a string.

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
Note that for the block syntax variables are not accepted because the evaluation of B is
done at parse time.
```
mybasis = 6-31g
@set {
    basis mybasis
}
```
Will set the basis to "mybasis" not "6-31g".
"""
macro set(opt,val)
    clean_up(s) = String(filter(c->!occursin(c," ():"),s))
    A = clean_up(repr(opt))
    B = clean_up(repr(val))
    quote
        try
            CurrentOptions[$A] = $val
        catch UndefVarError
            CurrentOptions[$A] = $B
        end
    end |> esc
end

macro set(block)
    lines = split(repr(block),";")
    quote 
        clean_up(s) = String(strip(filter(c->!occursin(c," {}():"),s)))
        for l in $lines
            opt, val = split(strip(l), " ", limit=2)
            opt = clean_up(opt)
            val = clean_up(val)
            try
                CurrentOptions[opt] = eval(Meta.parse(val))
            catch UndefVarError
                CurrentOptions[opt] = val
            end
        end
    end
end

macro get(opt)
    clean_up(s) = String(filter(c->!occursin(c," ():"),s))
    A = clean_up(repr(opt))
    quote
        try
            CurrentOptions[$A]
        catch KeyError
            throw(InvalidFermiOption($A*" not a valid option."))
        end
    end |> esc
end
"""
    Fermi.@molecule

Set the molecule used in computations. String is save into Fermi.CurrentOptions

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

struct InvalidFermiOption <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidFermiOption) = print(io, "InvalidFermiOption: ", e.msg)

#function notimplemented()
#    @output "🚧 Not implemented yet! We're working on it 🔨 👷 \n"
#end
