using Fermi
using Test

wf = @energy rhf;
@test wf.energy == -74.9650028737304410
@set basis cc-pvdz
@set scf_alg df
@energy rhf;
#@test wf.energy ≈
