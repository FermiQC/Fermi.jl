from pyscf import scf, gto

molstring = open("S22_15_dimer.xyz", "r").read()

mol = gto.M(atom = molstring[6:], basis="cc-pvdz", symmetry=False)

rhf = scf.RHF(mol)
rhf.direct_scf = False
rhf.density_fit(auxbasis="cc-pvtz-jkfit")
rhf.kernel()

print(rhf)
