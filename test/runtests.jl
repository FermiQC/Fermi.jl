using Fermi
using Test

wf = @energy rhf;
@test isapprox(wf.energy, -74.9650028737304410, rtol=1E-10)
wf = @energy ccsd;
@test isapprox(wf.CorrelationEnergy, -75.0187095714, rtol=1E-10)
@set basis cc-pvdz
@set scf_alg df
@energy rhf;
#@test wf.energy ≈
