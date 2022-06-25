using GaussianBasis
using Molecules
using TensorOperations

function RCCSDpTgrad(x...)
    RCCSDpTgrad(Molecule(), x...)
end

function RCCSDpTgrad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(FermiException("Invalid or unsupported derivative type for RCCSDpT: \"$dtype\""))
    elseif dtype == "findif"
        Fermi.gradient_findif(Fermi.CoupledCluster.RCCSDpT, mol, x...)
    else
        throw(FermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###