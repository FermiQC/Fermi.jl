# default options go here
# if user specifies option they will be overwritten
export CurrentOptions
export InvalidFermiOption
export @set

"""
    Fermi.CurrentOptions

Dictionary containing options for Fermi. Any information not given
explicitly to Methods is obtained from here.
"""
CurrentOptions = Dict{String,Any}(
                                  "molstring" => """
                                  O        1.2091536548      1.7664118189     -0.0171613972
                                  H        2.1984800075      1.7977100627      0.0121161719
                                  H        0.9197881882      2.4580185570      0.6297938832
                                  """,
                                  "basis" => "sto-3g",
                                  "charge" => 0,
                                  "multiplicity" => 1,
                                  "unit" => "angstrom",
                                  "reference" => "rhf",
                                  "scf_max_iter" => 50,
                                  "scf_max_rms" => 10^-10,
                                  "scf_alg" => "conventional",
                                  "quiet" => "true",
                                  "mp2_type" => "conv",
                                  "precision" => "double",
                                  "cc_alg" => "CTF",
                                  "ci_alg" => "sparse",
                                  "e_conv" => 10,
                                  "d_conv" => 8,
                                  "cc_max_iter" => 50,
                                  "cc_max_rms" => 10^-10,
                                  "cc_e_conv" => 10^-10,
                                  "drop_occ" => 0,
                                  "drop_vir" => 0,
                                  "diis" => false,
                                  "num_frozen" => 0,
                                  "cas_frozen" => 0,
                                  "cas_active" => -1,
                                  "cas_cutoff" => 10^-9,
                                  "cas_nroot" => 1,
                                  "min_matrix_elem" => 10^-9
                                 )

struct InvalidFermiOption <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidFermiOption) = print(io, "InvalidFermiOption: ", e.msg)


"""
    Fermi.@set

Set options for Fermi computations. It saves the options into Fermi.CurrentOptions.

*Usage:*  @set A B

A is set to B. By default A is taken as a string. B is evaluated at parse time. If the evaluation
is not possible, B is converted to a string.

# Examples

```
@set basis cc-pVDZ
@set cc_max_iter 100
@set e_conv 10^-9
```
"""
macro set(opt,val)
    clean_up(s) = filter(c->!occursin(c," ():"),s)
    A = clean_up(repr(opt))
    try
        CurrentOptions[A] = eval(val)
    catch UndefVarError
        CurrentOptions[A] = clean_up(repr(val))
    end
end
