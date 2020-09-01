export @energy

energy_dict = Dict{String, Expr}(
    "rhf" => :(Fermi.HartreeFock.RHF()),
    "rmp2" => :(Fermi.MollerPlesset.RMP2()),
    "rmp3" => :(Fermi.MollerPlesset.RMP3()),
    "ccsd" => :(Fermi.CoupledCluster.RCCSD()),
    "ecccsd" => :(Fermi.CoupledCluster.ecRCCSD()),
    "ccsdpt"=> :(Fermi.CoupledCluster.RCCSDpT()),
    "ecccsdpt"=> :(Fermi.CoupledCluster.ecRCCSDpT()),
    "ci" => :(Fermi.ConfigurationInteraction.CASCI()),
)

"""
    @energy(comm)

Macro to call functions to compute energy given current options.

# Implemented methods:
    rhf            Restricted Hartree-Fock
    rmp2           Restricted Moller-Plesset PT order 2
    rmp3           Restricted Moller-Plesset PT order 3
    ccsd           Restricted Coupled-Cluster with Single and Double substitutions
    ccsdpt         Restricted Coupled-Cluster with Single and Double substitutions and perturbative triples.
    ecCCSD         Restricted externally corrected CCSD. External correction from CASCI computation.
    ecCCSDpT       Restricted externally corrected CCSDpT. External correction from CASCI computation.
"""
macro energy(comm)
    clean_up(s) = lowercase(String(filter(c->!occursin(c," ():"),s)))
    A = clean_up(repr(comm))
    Fermi.energy_dict[A]
end
