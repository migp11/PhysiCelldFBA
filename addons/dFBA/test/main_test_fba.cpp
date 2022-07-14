#include <iostream>
#include <vector>

#include "../src/FBA_model.h"

using namespace std;

int main (int argc, const char *argv[])
{

    if (argc != 2)
    {
        cout << endl << "Usage: readSBML filename" << endl << endl;
        return 1;
    }
    const char *sbml_fileame = argv[1];

    FBA_model *model = new FBA_model();
    std::cout << "Reading SBML model from: " << sbml_fileame << " ";
    model->readSBMLModel(sbml_fileame);
    std::cout << "Ok!" << std::endl;
    std::cout << "Model " << model->getId() << "has been correctly loaded" << std::endl;
    std::cout << "Initializing LP model: ";
    model->initLpModel();
    std::cout << "Ok!" << std::endl;

    std::cout << "Testing FBA: ";
    FBA_solution* solution = model->optimize();
    if ( solution->status == 'optima'){
        std::cout << "OPTIMAL SOLUTION FOUND" << std::endl;
        FBA_solution* solution = model->getSolution();
        float fopt = solution->getObjectiveValue();
	    std::cout << "Objective value: " << fopt << std::endl;   
    }

    delete model;
    return 0;
}
