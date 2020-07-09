using Fermi

struct Molecule
    molstring::String
    charge::Int
    multiplicity::Int
    Nα::Int
    Nβ::Int
    Vnuc::Float64
end

function Molecule(mol::String, unit::String="angstrom")
    charge = Fermi.CurrentOptions["charge"]
    multiplicity = Fermi.CurrentOptions["multiplicity"]
    Molecule(mol, charge, multiplicity, unit)
end

function Molecule(mol::String, charge::Int, multiplicity::Int, unit::String="angstrom")
    
    if unit == "bohr"
        Vnuc = nuclear_repulsion(mol, unit="bohr")
        molstring = convert_xyz_unit(mol, from="bohr", to="angstrom")
    elseif unit == "angstrom"
        Vnuc = nuclear_repulsion(mol, unit="angstrom")
        molstring = mol
    end
    
    Nα, Nβ = get_num_electrons(molstring, charge, multiplicity)

    return Molecule(molstring, charge, multiplicity, Nα, Nβ, Vnuc)
end

function convert_xyz_unit(molstring::String; from::String, to::String)

    if from == "angstrom" && to == "bohr"
        conv = Fermi.PhysicalConstants.angstrom_to_bohr
    elseif from == "bohr" && to == "angstrom"
        conv = Fermi.PhysicalConstants.bohr_to_angstrom
    else
        error("Invalid unit convertion requested at Fermi.Integrals.convert_xyz_unit. Requested: $from → $to. Supported: Bohr ↔ Angstrom")
    end

    new_molstring = ""
    exp = r"(\w{1,2})\s+([-+]??\d+\.\d+\s+[-+]??\d+\.\d+\s+[-+]??\d+\.\d+)"
    for line in split(molstring, "\n")
        m = match(exp, line)
        if m != nothing

            # Convert SubString to String
            A = String(m.captures[1])

            # Convert Substring to String and then String to Float64 Array
            P = parse.(Float64, split(String(m.captures[2]))).*conv

            # Save atomic numbers and positions
            new_molstring *= A
            new_molstring *= "   "
            new_molstring *= "$(P[1])   $(P[2])   $(P[3])\n"
        end
    end
    return String(strip(new_molstring))
end

function nuclear_repulsion(mol::String, unit::String = "angstrom")

    # Check unit
    if unit == "angstrom"
        molstring = convert_xyz_unit(mol, from="angstrom", to="bohr")
    elseif unit == "bohr"
        molstring = mol
    else
        throw(Fermi.InvalidFermiOption("$unit is not a known keyword for unit"))
    end

    Zvals = Int64[]
    positions = Array{Float64,1}[]

    # Reg exp for xyz lines
    exp = r"(\w{1,2})\s+([-+]??\d+\.\d+\s+[-+]??\d+\.\d+\s+[-+]??\d+\.\d+)"
    for line in split(molstring, "\n")
        m = match(exp, line)
        if m != nothing

            # Convert SubString to String
            A = String(m.captures[1])

            # Convert Substring to String and then String to Float64 Array
            P = parse.(Float64, split(String(m.captures[2])))

            # Save atomic numbers and positions
            push!(Zvals, Fermi.PhysicalConstants.atomic_number(A))
            push!(positions, P)
        end
    end

    # Compute Nuclear repulsion
    Vnuc = 0.0
    for i in eachindex(Zvals)
        @inbounds Zi = Zvals[i]
        @inbounds A = positions[i]
        for j in 1:(i-1) # Avoid double counting
            @inbounds Zj = Zvals[j]
            @inbounds B = positions[j]
            Vnuc += (Zi*Zj)/√((A.-B)⋅(A.-B))
        end
    end
    return Vnuc
end

function get_num_electrons(molstring::String, charge::Int, multiplicity::Int)

    nelec = 0
    # Reg exp for xyz lines
    exp = r"(\w{1,2})\s+([-+]??\d+\.\d+\s+[-+]??\d+\.\d+\s+[-+]??\d+\.\d+)"
    for line in split(molstring, "\n")
        m = match(exp, line)
        if m != nothing

            # Convert SubString to String
            A = String(m.captures[1])

            # Save atomic numbers and positions
            nelec += Fermi.PhysicalConstants.atomic_number(A)
        end
    end

    nelec -= charge
    αexcess = multiplicity-1

    if isodd(nelec) != isodd(αexcess)
        throw(Fermi.InvalidFermiOption("Incompatible charge $(charge) and multiplicity $(multiplicity)"))
    end

    Nβ = (nelec - αexcess)/2
    Nα = nelec - Nβ

    return Nα, Nβ
end
