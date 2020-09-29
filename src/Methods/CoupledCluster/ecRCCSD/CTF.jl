"""
    Fermi.CoupledCluster.RCCSD{T}(Alg::CTF)

Compute a RCCSD wave function using the Compiled time factorization algorithm (CTF)
"""
function ecRCCSD{T}(Alg::CTF) where T <: AbstractFloat
    @output "Calling CASCI module...\n"
    # Call CASCI
    cas = Fermi.ConfigurationInteraction.CASCI()
    ecRCCSD{T}(cas, Alg)
end

function ecRCCSD{T}(refwfn::Fermi.HartreeFock.RHF, Alg::CTF) where T <: AbstractFloat
    @output "Using given RHF to compute CASCI...\n"
    cas = Fermi.ConfigurationInteraction.CASCI(refwfn)
    ecRCCSD{T}(cas, Alg)
end

function ecRCCSD{T}(refwfn::Fermi.HartreeFock.RHF, cas::Fermi.ConfigurationInteraction.CASCI, Alg::CTF) where T <: AbstractFloat

    # Save reference wavefunction and process CAS data. Modify Ref (if not HF)
    @output "\n\nProcessing CAS data...\n"
    refdet, Casdata = process_cas(cas)
    ecRCCSD{T}(refdet, refwfn, Casdata, Alg)
end

function ecRCCSD{T}(cas::Fermi.ConfigurationInteraction.CASCI, Alg::CTF) where T <: AbstractFloat

    # Save reference wavefunction and process CAS data. Modify Ref (if not HF)
    @output "\n\nProcessing CAS data...\n"
    refdet, Casdata = process_cas(cas)
    refwfn = cas.ref
    ecRCCSD{T}(refdet, refwfn, Casdata, Alg)
end

function ecRCCSD{T}(refdet::Determinant, refwfn::Fermi.HartreeFock.RHF, Casdata::Dict{String,Array}, Alg::CTF) where T <: AbstractFloat
    # Print intro
#    Fermi.CoupledCluster.print_header()
    @output "\n    • Computing Externally Corrected CCSD with the ecRCCSD module.\n\n"

    # Recover integrals object
    ints = refwfn.ints

    # Get MO Integrals
    drop_occ = Fermi.CurrentOptions["drop_occ"]
    drop_vir = Fermi.CurrentOptions["drop_vir"]

    # Process CAS data to get ecT1 and ecT2 (Cluster Decomposition step)
    frozen = Fermi.CurrentOptions["cas_frozen"]
    active = Fermi.CurrentOptions["cas_active"] ≢ -1 ? Fermi.CurrentOptions["cas_active"] : refwfn.nvir+refwfn.ndocc-frozen

    if drop_occ > frozen
        error("\nFrozen orbitals in the CC step ($drop_occ) cannot be greater than the number of frozen orbitals in the CAS ($frozen).")
    end

    if drop_vir > (refwfn.ndocc+refwfn.nvir) - active - frozen
        error("\nToo many virtual orbitals dropped ($drop_vir) for the active space.")
    end

    @output "\n   • CAS Decomposition started:\n"
    t = @elapsed T1, T2, ecT1, ecT2 = cas_decomposition(refdet, Casdata, refwfn.ndocc, drop_occ, ints["FOV"], ints["OOVV"], ints["OVVV"], ints["OOOV"])
    @output "Finished in {:5.5} seconds.\n" t
    Casdata = nothing

    RCCSD{T}(refwfn, ints, T1, T2, Alg, ecT1=ecT1, ecT2=ecT2)
end

function process_cas(cas::Fermi.ConfigurationInteraction.CASCI)
    
    # This function process the CAS wave function to return a dictionary with C coefficients 
    # organized by excitation level, along with the corresponding determinant
    # CI coefficients are intermediate normalized.
    # The reference determinant is taken as the HF one.

    dets = cas.dets
    Ccas = cas.coef

    ref = Determinant(repeat("1",cas.ref.ndocc),repeat("1",cas.ref.ndocc))
    z = 0
    acf = 0

    @output "   • CAS Composition\n"
    @output "Excitation      N of dets\n"
    while acf < length(dets)
        x = count(d->excitation_level(ref,d)==z, dets)
        @output "{}           {}\n" z x
        z += 1
        acf += x
    end

    ref = dets[1]
    zeroth = repeat('1', cas.ref.ndocc)*repeat('0', cas.ref.nvir)
    hf = Fermi.ConfigurationInteraction.Determinant(zeroth, zeroth)
    if ref == hf
        @output "Dominant configuration is the RHF determinant.\n"
    else
        error("Dominant determinant is not the RHF.")
    end
    C0 = Ccas[1]

    # Intermediate Normalization
    Ccas = Ccas ./ C0

    # Split the Cas data into excitation level
    Casdata = Dict{String,Array}()

    Casdata["C1"] = Float64[]
    Casdata["C2"] = Float64[]
    Casdata["C3"] = Float64[]
    Casdata["C4"] = Float64[]

    Casdata["C1dets"] = Determinant[]
    Casdata["C2dets"] = Determinant[]
    Casdata["C3dets"] = Determinant[]
    Casdata["C4dets"] = Determinant[]

    for i in eachindex(dets)

        exc = excitation_level(ref, dets[i])

        if exc == 1 
            push!(Casdata["C1"], Ccas[i])
            push!(Casdata["C1dets"], dets[i])

        elseif exc == 2
            push!(Casdata["C2"], Ccas[i])
            push!(Casdata["C2dets"], dets[i])

        elseif  exc == 3
            push!(Casdata["C3"], Ccas[i])
            push!(Casdata["C3dets"], dets[i])

        elseif  exc == 4
            push!(Casdata["C4"], Ccas[i])
            push!(Casdata["C4dets"], dets[i])
        end
    end
    return ref, Casdata
end
