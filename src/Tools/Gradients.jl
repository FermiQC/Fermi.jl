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
using "=>" or "<=" for analytic gradients. The derivative type (analytic, findif, autodif) 
can be set with `@set deriv_type <kw>`.

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