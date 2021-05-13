#include <iostream>
#include <vector>

#include "../src/FBA_model.h"

#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>

using namespace std;
using namespace FBA;

int main(int argc, const char *argv[])
{

    if (argc != 2)
    {
        cout << endl << "Usage: readSBML filename" << endl << endl;
        return 1;
    }
    const char* sbml_fileame   = argv[1];

    FBA_model* model = new FBA_model();
    std::cout << "Reading SBML model from: " << sbml_fileame;
    model->readSBMLModel(sbml_fileame);
    std::cout << "Ok!" << std::endl;
    std::cout << "Reading SBML model from: " << sbml_fileame << std::endl;

    std::vector<FBA_reaction*> listOfReactions = model->getListOfReactions();
    for(FBA_reaction* r: listOfReactions) {
        std::string rId = r->getId();
        double lb = r->getLowerBound();
        double ub = r->getUpperBound();
        double c = r->getObjectiveCoefficient();
        std::cout << rId << " c=" << c << " (" << lb << "," << ub << ")";
        std::cout << std::endl;
        std::cout << "  " << r->getReactionString() << std::endl;
        double Km = r->getKm();
        double Vmax = r->getVmax();
        double substrate = 10;
        double bound = r->calculateFluxBound(substrate);
        std::cout << "Km: " << Km << " ";
        std::cout << "Vmax: " << Vmax << " ";
        std::cout << "Bound: " << bound << "\n\n";
    }

    delete model;
    return 0;
}
