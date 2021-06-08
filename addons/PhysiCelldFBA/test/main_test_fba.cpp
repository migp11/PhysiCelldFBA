#include <iostream>
#include <vector>

#include "../src/FBA_model.h"

using namespace std;

#define ZERO_TOLERANCE 1e-7;

void printSolution(FBA_model* model)
{
    for( FBA_reaction* exchange: model->getListOfBoundaryReactions() )
    {
        float flux = exchange->getFluxValue();
        if( abs(flux) < 1e-7 )
            continue;
            
        std::cout << exchange->getId() << ": " << flux << endl;
    }
}

int
main (int argc, const char *argv[])
{

    
    string exchange_o2 = "R_EX_o2_e";
    string exchange_glc = "R_EX_glc__D_e";

    const char *sbml_fileame;
    if( argc > 1 )
	    sbml_fileame = argv[1];
	else
	    sbml_fileame = "test/data/Ecoli_core.xml";

    cout << "=========================================" << endl;
    cout << " TEST FBA_models package (libSBML + Clp) " << endl;
    cout << "=========================================" << endl;
    
    FBA_model model;
    
    std::cout << "Loading SBML model from: " << sbml_fileame << " ";
    model.readSBMLModel(sbml_fileame);
    std::cout << "Ok!" << std::endl;
    
    std::cout << "Initializing LP model: ";
    model.initLpModel();
    std::cout << "Ok!" << std::endl;

    cout << endl;
    cout << "=============================" << endl;
    cout << " Running first test (aerobic) " << endl;
    cout << "=============================" << endl;
        
    model.runFBA();
    std::cout << "Running FBA: ";
    cout << endl;
    if (model.getSolutionStatus ())
    {
        std::cout << "OPTIMAL SOLUTION FOUND" << std::endl;
        cout << endl;
        float fopt = model.getObjectiveValue();
        std::cout << "Objective value: " << fopt << std::endl;
        printSolution(&model);
    }
    else
    {
        std::cout << "NO SOLUTION FOUND" << std::endl;
    }
    
    cout << endl;
    cout << "================================" << endl;
    cout << " Running second test (anaerobic) " << endl;
    cout << "================================" << endl;
    
    float lb = model.getReactionLowerBound(exchange_o2);
    model.setReactionLowerBound(exchange_o2, 0.0);
    std::cout << "- Changing " << exchange_o2 << " lb from " << lb << " to " << 0.0 << std::endl;

    model.runFBA();
    std::cout << "Running FBA: ";
    if (model.getSolutionStatus ())
    {
        std::cout << "OPTIMAL SOLUTION FOUND" << std::endl;
        cout << endl;
        float fopt = model.getObjectiveValue();
        std::cout << "Objective value: " << fopt << std::endl;
        printSolution(&model);
    }
    else
    {
        std::cout << "NO SOLUTION FOUND" << std::endl;
    }
    

    cout << endl;
    cout << "======================================" << endl;
    cout << " Running third test test (infeasible) " << endl;
    cout << "======================================" << endl;
    
    lb = model.getReactionLowerBound(exchange_glc);
    model.setReactionLowerBound(exchange_glc, 0.0);
    std::cout << "- Changing " << exchange_glc << " lb from " << lb << " to " << 0.0 << std::endl;

    model.runFBA();
    std::cout << "Running FBA: ";
    if (model.getSolutionStatus ())
    {
        std::cout << "OPTIMAL SOLUTION FOUND" << std::endl;
        cout << endl;
        float fopt = model.getObjectiveValue();
        std::cout << "Objective value: " << fopt << std::endl;
        printSolution(&model);
    }
    else
    {
        std::cout << "NO SOLUTION FOUND" << std::endl;
        printSolution(&model);
    }

    cout << endl;
    cout << "=============================" << endl;
    cout << " Test completed successfully! " << endl;
    cout << "=============================" << endl;
    
    return 0;
}
