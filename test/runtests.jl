using Fermi
using Test

wf = @energy rhf;
@test wf.energy == -74.9650028737304410
