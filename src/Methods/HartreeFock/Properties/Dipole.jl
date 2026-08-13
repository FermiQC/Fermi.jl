using GaussianBasis
using Molecules

# 1 a.u. of electric dipole moment (e·a₀) in Debye.
const AU2DEBYE = 2.541746

"""
    Fermi.HartreeFock.dipole(wfn::RHF)
    Fermi.HartreeFock.dipole(wfn::UHF)

Compute the electric dipole moment (electronic + nuclear) of a converged
Hartree-Fock wave function. Returns a 3-element vector `[μx, μy, μz]` in
atomic units (e·a₀).
"""
function dipole(wfn::RHF)
    atoms = wfn.molecule.atoms
    bset = BasisSet(wfn.orbitals.basis, atoms)

    # Dipole integrals from GaussianBasis are already in atomic units (Bohr),
    # since libcint is fed Bohr-converted coordinates internally. Molecules'
    # nuclear_dipole works directly off Atom.xyz, which is stored in Angstrom,
    # so it needs the same Å->Bohr conversion before it can be added to the
    # electronic term.
    dip = GaussianBasis.dipole(bset)

    Co = @view wfn.orbitals.C[:, 1:wfn.ndocc]
    D = Co*Co' # Not doubled -- see RHFEnergy/RHFgrad.jl for the same convention

    # Electron charge is -1 in atomic units, and D here counts one electron per
    # doubly-occupied spatial orbital, hence the extra factor of 2 (unlike the
    # RHF energy expression, this factor does not cancel against anything).
    μ_elec = [-2*sum(D .* view(dip, :, :, k)) for k = 1:3]
    μ_nuc = Molecules.nuclear_dipole(wfn.molecule) ./ Molecules.bohr_to_angstrom
    μ = μ_elec .+ μ_nuc

    print_dipole(μ)
    return μ
end

function dipole(wfn::UHF)
    atoms = wfn.molecule.atoms
    bset = BasisSet(wfn.orbitals.basis, atoms)
    dip = GaussianBasis.dipole(bset)

    Nα = wfn.molecule.Nα
    Nβ = wfn.molecule.Nβ
    Coα = @view wfn.orbitals.Cα[:, 1:Nα]
    Coβ = @view wfn.orbitals.Cβ[:, 1:Nβ]
    Dtot = Coα*Coα' .+ Coβ*Coβ' # Each spin already singly-occupied, no extra factor

    μ_elec = [-sum(Dtot .* view(dip, :, :, k)) for k = 1:3]
    μ_nuc = Molecules.nuclear_dipole(wfn.molecule) ./ Molecules.bohr_to_angstrom
    μ = μ_elec .+ μ_nuc

    print_dipole(μ)
    return μ
end

function print_dipole(μ)
    mag = √sum(abs2, μ)
    output("\n   • Dipole Moment")
    output(" {:>10} {:>15.10f} {:>15.10f}  (a.u., Debye)", "X", μ[1], μ[1]*AU2DEBYE)
    output(" {:>10} {:>15.10f} {:>15.10f}  (a.u., Debye)", "Y", μ[2], μ[2]*AU2DEBYE)
    output(" {:>10} {:>15.10f} {:>15.10f}  (a.u., Debye)", "Z", μ[3], μ[3]*AU2DEBYE)
    output(" {:>10} {:>15.10f} {:>15.10f}  (a.u., Debye)", "Total", mag, mag*AU2DEBYE)
end
