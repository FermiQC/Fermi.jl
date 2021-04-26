using Fermi
using Fermi.Error

const LIBPATH = joinpath(@__DIR__, "../../deps/lib")
const AM_pat = r"([SPDFGHI]{1,2})\s+?(\d++)"
const prim_pat = r"([+-]?\d*?\.\d+[D+-]{0,2}\d\d)\s+?([+-]?\d*?\.\d+[D+-]{0,2}+\d\d)"
const prim_pat3 = r"([+-]?\d*?\.\d+[D+-]{0,2}\d\d)\s+?([+-]?\d*?\.\d+[D+-]{0,2}+\d\d)\s+?([+-]?\d*?\.\d+[D+-]{0,2}+\d\d)"
const AMDict = Dict(
        "S" => 0,
        "P" => 1, 
        "D" => 2, 
        "F" => 3, 
        "G" => 4, 
        "H" => 5, 
        "I" => 6, 
    )

function read_basisset(bname::String, AtomSymbol::String)

    clean_bname = replace(bname, "*"=>"_st_")
    file_path = joinpath(LIBPATH, clean_bname*".gbs")

    output("\nLooking for atom $AtomSymbol in basisset file at: $file_path")

    if !(isfile(file_path))
        throw(InvalidFermiOption("Basis set file for $bname was not found."))
    end

    info = extract_atom_from_bs(file_path, AtomSymbol)

    # Split info into basis
    BasisStrings = [info[1]*"\n"]
    for line in info[2:end]
        if occursin(AM_pat, line)
            push!(BasisStrings, line*"\n")
        else
            BasisStrings[end] *= line * "\n"
        end
    end

    out = BasisFunction[]
    output("Basis functions for $AtomSymbol: ", ending="")
    for b in BasisStrings
        r = r"[SPDFGHI]{2}"
        if occursin(r, b)
            b1, b2 = two_basis_from_string(b)
            push!(out, b1)
            push!(out, b2)
        else
            push!(out, basis_from_string(b))
        end
    end
    output("")

    return out
end

function extract_atom_from_bs(file_path::String, AtomSymbol::String)

    atom_pat = Regex("$(AtomSymbol)\\s+?0")
    flag_atom = false

    out = String[]
    for line in eachline(file_path)

        if occursin(atom_pat, line)
            flag_atom = true
            continue
        end

        if flag_atom
            occursin(raw"****", line) ? break : nothing
            push!(out, String(line))
        end
    end

    if flag_atom
        return out
    else
        output("")
        throw(InvalidFermiOption("Atom $AtomSymbol not found in $file_path."))
    end
end

function basis_from_string(bstring::String)
    lines = split(strip(bstring), "\n")
    head = lines[1]
    m = match(AM_pat, head)

    if m === nothing
        throw(BasisSetError("cannot parse basis set file line: $head"))
    end

    AMsymbol = String(m.captures[1])
    l = AMDict[AMsymbol]
    nprim = parse(Int, m.captures[2])
    coef = zeros(nprim)
    exp = zeros(nprim)

    for i in eachindex(exp)
        m = match(prim_pat, lines[i+1])

        if m === nothing
            throw(BasisSetError("cannot parse basis set file line: $(lines[i+1])"))
        end

        e = replace(m.captures[1], "D"=>"E")
        c = replace(m.captures[2], "D"=>"E")
        coef[i] = parse(Float64, c)
        exp[i] = parse(Float64, e)
    end

    output(AMsymbol, ending=" ")
    return BasisFunction(l, coef, exp)
end

function two_basis_from_string(bstring::String)
    lines = split(strip(bstring), "\n")
    head = lines[1]
    m = match(AM_pat, head)
    if m === nothing
        throw(BasisSetError("cannot parse basis set file line: $head"))
    end
    AMsymbol = String(m.captures[1])

    length(AMsymbol) == 2 || throw(MethodArgument("cannot extract two basis from $AMsymbol function"))

    l1 = AMDict[AMsymbol[1]*""]
    l2 = AMDict[AMsymbol[2]*""]

    nprim = parse(Int, m.captures[2])
    coef1 = zeros(nprim)
    coef2 = zeros(nprim)
    exp = zeros(nprim)

    for i in eachindex(exp)
        m = match(prim_pat3, lines[i+1])

        if m === nothing
            throw(BasisSetError("cannot parse basis set file line: $(lines[i+1])"))
        end

        e = replace(m.captures[1], "D"=>"E")
        c1 = replace(m.captures[2], "D"=>"E")
        c2 = replace(m.captures[3], "D"=>"E")
        coef1[i] = parse(Float64, c1)
        coef2[i] = parse(Float64, c2)
        exp[i] = parse(Float64, e)
    end

    output(AMsymbol, ending=" ")
    return BasisFunction(l1, coef1, exp), BasisFunction(l2, coef2, exp)
end