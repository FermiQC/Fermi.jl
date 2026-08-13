export @gradient

g_dict = Dict{String, String}(
    "rhf" => "Fermi.HartreeFock.RHFgrad()",
    "uhf" => "Fermi.HartreeFock.UHFgrad()",
    "mp2" => "Fermi.MollerPlesset.RMP2grad()",
    "rmp2" => "Fermi.MollerPlesset.RMP2grad()",
    "ccsd" => "Fermi.CoupledCluster.RCCSDgrad()",
    "rccsd" => "Fermi.CoupledCluster.RCCSDgrad()",
    "ccsd(t)"=> "Fermi.CoupledCluster.RCCSDpTgrad()",
    "rccsd(t)"=> "Fermi.CoupledCluster.RCCSDpTgrad()"
)

"""
    Fermi.@gradient

Macro to call functions to compute gradients given current options. Arguments may be passed
using "=>" or "<=" for analytic gradients. `@set deriv_type <kw>` controls how the gradient
is obtained:

    "auto"      (default) Analytic if the method implements one (currently only RHF),
                otherwise finite differences. Never errors due to missing analytic support.
    "analytic"  Force analytic. Throws if the method has no analytic gradient.
    "findif"    Force finite differences, even for methods with an analytic gradient
                available (e.g. to cross-check RHF's analytic gradient).

# Examples

Generating a RHF gradient
```
@gradient rhf
```

By default, gradients are calculated using the current `molstring`, but `Molecule` objects
can also be passed to the gradients
```
mol = Molecule(molstring=mymol)
@gradient mol => rhf
```
or
```
@gradient rhf <= mol
```

A bare `@gradient rhf` (or `mol => rhf`) has no wave function or atomic
integrals (`aoints`) to reuse, so it runs SCF from scratch and prints a
warning before doing so (real, unavoidable work -- SCF needs its atomic
integrals regardless). To avoid that when you already have a wave function
at the same geometry (e.g. right after `@energy`), pass it in directly --
`@gradient wfn => rhf` builds its own throwaway `aoints` internally, but
prints no warning: the gradient itself never touches `aoints`'s cache, only
cheap basis/molecule metadata, so nothing expensive is actually being
redone. That internal `aoints` still follows the current `@set df <bool>`
option though (same as a fresh `IntegralHelper()` would), so the gradient
matches whatever ERI representation `wfn` was actually converged with.

For full control (e.g. mixing `df` settings between energy and gradient, or
truly avoiding any extra `IntegralHelper` construction), pass `aoints`
explicitly:
```
aoints = Fermi.Integrals.IntegralHelper()
wfn = @energy aoints => rhf
grad = @gradient aoints, wfn => rhf
```
RHF's analytic gradient follows whatever `eri_type` `aoints` carries --
`Chonky`/`SparseERI` use the exact ERI, while a density-fitted `aoints`
(`@set df true`, or an explicit `JKFIT()`/`RIFIT()` `eri_type`) gets the
matching RI-HF analytic gradient (no CPHF needed, same as the exact-ERI
case).

# Implemented methods:
    Method   Type      Description
    ------   ----      -----------
    RHF      AN,FD     Restricted Hartree-Fock.
    UHF      FD        Unrestricted Hartree-Fock
    RMP2     FD        Restricted Moller-Plesset PT order 2.
    CCSD     FD        Restricted Coupled-Cluster with Single and Double substitutions
    CCSD(T)  FD        Restricted Coupled-Cluster with Single and Double substitutions and perturbative triples.
"""
macro gradient(comm)
    # Clean spaces and colon from string
    clean_up(s) = String(filter(c->!occursin(c," :"),s))

    # Transform the input expression into a string
    A = clean_up(repr(comm))

    # If the command string is in between parenthesis, remove it
    while A[1] == '(' && A[end] == ')'
        A = A[2:end-1]
    end

    # If "=>" is present, then the second part is taken as the method and the first as the argument passed
    if occursin("=>", A)
        arg, method = split(A, "=>")
    elseif occursin("<=", A)
        method, arg = split(A, "<=")
    else
        method = A
        arg = ""
    end

    method = lowercase(String(method))
    arg = String(arg)   # Make sure arg is not a substring
    
    out = ""
    try
        out = replace(g_dict[method], "()" => "("*arg*")")
    catch KeyError
        throw(FermiException("Invalid or unsupported method for gradient computation: \"$A\""))
    end

    expr_out = Meta.parse(out)

    return esc(expr_out)
end