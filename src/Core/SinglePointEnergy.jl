export @energy

energy_dict = Dict{String, Expr}(
    "rhf" => :(Fermi.HartreeFock.RHF()),
    "ccsd" => :(Fermi.CoupledCluster.RCCSD()),
    "ccsd(t)"=> :(Fermi.CoupledCluster.RCCSDpT()),
    "ci" => :(Fermi.ConfigurationInteraction.CASCI()),
)

macro energy(comm)
    clean_up(s) = lowercase(String(filter(c->!occursin(c," ():"),s)))
    A = clean_up(repr(comm))
    Fermi.energy_dict[A]
end
