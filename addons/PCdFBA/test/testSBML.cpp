#include <iostream>
#include <string>


/*
comment from post #define LIBSBML_STATIC #define LIBLAX_STATIC
*/

#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>

LIBSBML_CPP_NAMESPACE_USE


using namespace std;

int main (int argc, char* argv[])
{

  SBMLReader reader;
  SBMLDocument* document;

  const char *sbml_fileame;
  if( argc > 1 )
	    sbml_fileame = argv[1];
	else
	    sbml_fileame = "test/data/Ecoli_core.xml";

  cout << endl;
  cout << "=================" << endl;
  cout << " TEST libSBML" << endl;
  cout << "=================" << endl;
  
  document = reader.readSBML(sbml_fileame);
  unsigned int errors = document->getNumErrors();

  cout << endl;
  cout << " Reading filename: " << sbml_fileame << endl;
  cout << " Validation error(s): " << errors << endl;
  cout << endl;
  document->printErrors(cerr);

  Model* model = document->getModel();
  
  int n_species = model->getNumSpecies();
  int n_reactions = model->getNumReactions();
  int n_parmas = model->getNumParameters();
  string model_id = model->getId();
  
  cout << " Model ID: " << model_id << endl;
  cout << " - Total species: " << n_species << endl;
  cout << " - Total reactions: " << n_reactions << endl;
  cout << " - Total params: " << n_parmas << endl;

  cout << endl;
  cout << "=============================" << endl;
  cout << " Test completed successfully! " << endl;
  cout << "=============================" << endl;

  delete document;
  
  return errors;
}


