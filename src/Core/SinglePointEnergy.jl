export @energy

energy_dict = Dict{String, String}(
    "rhf" => "Fermi.HartreeFock.RHF()",
)

"""
    Fermi.@energy

Macro to call functions to compute energy given current options. Arguments may be passed
using "=>" or "<="

# Examples

Generating a RHF wave function
```
rhf = @energy rhf
```

Using this wave function in a CCSD computations
```
@energy rhf => CCSD
```
or
```
@energy CCSD <= rhf
```

# Implemented methods:
    RHF            Restricted Hartree-Fock.
    RMP2           Restricted Moller-Plesset PT order 2.
    RMP3           Restricted Moller-Plesset PT order 3.
    CASCI          Restricted configuration interaction within given active space.
    CCSD           Restricted Coupled-Cluster with Single and Double substitutions
    CCSD(T)        Restricted Coupled-Cluster with Single and Double substitutions and perturbative triples.
    ecCCSD         Restricted externally corrected CCSD. External correction from CASCI computation.
    ecCCSD(T)      Restricted externally corrected CCSDpT. External correction from CASCI computation.
"""
macro energy(comm)
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
         out = replace(energy_dict[method], "()" => "("*arg*")")
    catch KeyError
        throw(Fermi.InvalidFermiOption("Invalid method for energy computation: \"$A\""))
    end

    expr_out = Meta.parse(out)

    return esc(expr_out)
end
