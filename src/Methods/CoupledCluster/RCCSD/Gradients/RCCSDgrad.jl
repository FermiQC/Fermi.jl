using GaussianBasis
using Molecules
using TensorOperations

function RCCSDgrad(x...)
    RCCSDgrad(Molecule(), x...)
end

function RCCSDgrad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(FermiException("Invalid or unsupported derivative type for RCCSD: \"$dtype\""))
    elseif dtype == "findif"
        Fermi.gradient_findif(Fermi.CoupledCluster.RCCSD, mol, x...)
    else
        throw(FermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###