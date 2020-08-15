using Fermi
using Test

tol = 1E-8
# RHF 
wf = @energy rhf;
@test isapprox(wf.energy, -74.9650028737304410, rtol=tol)
# CCSD
wf = @energy ccsd;
@test isapprox(wf.CorrelationEnergy, -75.0187095714, rtol=tol)
# externally corrected CC
wf = @energy ecCCSD;
@test isapprox(wf.CorrelationEnergy, -75.0187251416, rtol=tol)
#@set basis cc-pvdz
#@set scf_alg df
#@energy rhf;
#@molecule {
#          O 0.0 0.0 0.0
#          H 1.0 0.0 0.0
#          H 0.0 1.0 0.0
#         }
#@set cc_alg CTF
#wf = @energy ecCCSD
#@test isapprox(wf.CorrelationEnergy, -76.2329543810, rtol=tol)
