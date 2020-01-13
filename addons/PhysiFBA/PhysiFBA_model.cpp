/*
 * flux_balance_model.cpp
 *
 *  Created on: 11 jun. 2019
 *      Author: mponce
 */

#include "PhysiFBA_model.h"


namespace PhysiFBA{

  /* Default model used to initialize the initial cell */
  PhysiFBA_model PhysiFBA_default_model;
  std::map<std::string, std::string> exchange_flux_density_map;


  PhysiFBA_model::PhysiFBA_model()
  {
        this->id = "None";
        this->lp_model = nullptr;
  }

  PhysiFBA_model::~PhysiFBA_model() {
    for(PhysiFBA_reaction* rxn: this->reactions)
         delete rxn;

     for(PhysiFBA_metabolite* met: this->metabolites)
         delete met;

     delete lp_model;

  }

  const ClpSimplex* PhysiFBA_model::getLpModel() const
  {
      return this->lp_model;
  }

  const int PhysiFBA_model::getNumReactions()
  {

      return this->reactions.size();
  }

  const int PhysiFBA_model::getNumMetabolites()
  {
      return this->metabolites.size();
  }

  bool PhysiFBA_model::hasMetabolite(std::string mId)
  {
      std::map<std::string, int>::iterator itr;
      itr = this->metaboliteIndexer.find(mId);
      return itr != this->metaboliteIndexer.end();
  }

  void PhysiFBA_model::addMetabolite(PhysiFBA_metabolite* met)
  {
      if (!this->hasMetabolite( met->getId() ))
      {
          this->metabolites.push_back(met);
          this->metaboliteIndexer[met->getId()] = this->metabolites.size() - 1;
      }

  }

  const PhysiFBA_metabolite* PhysiFBA_model::getMetabolite(std::string mId)
  {
      if (this->hasMetabolite(mId))
      {
          int idx = this->metaboliteIndexer[mId];
          return this->metabolites[idx];
      }
      return nullptr;
  }

  const std::vector<PhysiFBA_metabolite*> PhysiFBA_model::getListOfMetabolites() const
  {
      return this->metabolites;
  }

  bool PhysiFBA_model::hasReaction(std::string rId)
  {
      std::map<std::string, int>::iterator it;
      it = this->reactionsIndexer.find(rId);

      return it != this->reactionsIndexer.end();
  }

  PhysiFBA_reaction* PhysiFBA_model::getReaction(std::string rId)
  {
      if (this->hasReaction(rId))
      {
          int idx = this->reactionsIndexer[rId];
          PhysiFBA_reaction* rxn = this->reactions[idx];
          return rxn;
      }
      return nullptr;
  }

  void PhysiFBA_model::setReactionUpperBound(std::string rId, float upperBound)
  {
    PhysiFBA_reaction* rxn = this->getReaction(rId);
    if (rxn)
      {
        rxn->setUpperBound(upperBound);
        int colIdx = this->reactionsIndexer[rId];
        this->lp_model->setColumnUpper(colIdx, upperBound);
      }

  }
  void PhysiFBA_model::setReactionLowerBound(std::string rId, float lowerBound)
  {
    PhysiFBA_reaction* rxn = this->getReaction(rId);
    if (rxn)
      {
        rxn->setLowerBound(lowerBound);
        int colIdx = this->reactionsIndexer[rId];
        this->lp_model->setColumnLower(colIdx, lowerBound);
      }

  }

  void PhysiFBA_model::addReaction(PhysiFBA_reaction* rxn)
  {
      if (!this->hasReaction( rxn->getId() ))
      {
          this->reactions.push_back(rxn);
          this->reactionsIndexer[rxn->getId()] = this->reactions.size() - 1;
      }
  }

  const int PhysiFBA_model::getReactionIndex(std::string rId)
  {
      if (this->hasReaction(rId))
          return this->reactionsIndexer[rId];
      else
          return -1;
  }


  const std::vector<PhysiFBA_reaction*> PhysiFBA_model::getListOfReactions() const
  {
      return this->reactions;
  }

  std::vector<std::string> PhysiFBA_model::getListOfBoundaryReactionIds()
  {
      std::vector<std::string> listOfBoundaryIds;
      for(PhysiFBA_reaction* reaction: this->reactions)
      {
          if (reaction->getNumberOfMetabolites() == 1)
          {
              listOfBoundaryIds.push_back(reaction->getId());
          }
      }
      return listOfBoundaryIds;
  }

  void PhysiFBA_model::readSBMLModel(const char* sbmlFileName)
  {
      libsbml::SBMLReader reader;
      libsbml::SBMLDocument* document = reader.readSBML(sbmlFileName);
      libsbml::Model* model = document->getModel();

      libsbml::ListOfSpecies* listOfSpecies = model->getListOfSpecies();
      libsbml::ListOfReactions* listOfPhysiFBA_reactions = model->getListOfReactions();
      libsbml::ListOfParameters* listOfParameters = model->getListOfParameters();

      this->id = model->getId();

      for (unsigned int i = 0; i < model->getNumSpecies(); i++)
      {
          libsbml::Species* species = listOfSpecies->get(i);
          // Skipping boundary metabolites
          if ( species->getBoundaryCondition() )
              continue;

          PhysiFBA_metabolite* metabolite = new PhysiFBA_metabolite(species->getId());
          metabolite->setName(species->getName());
          this->addMetabolite(metabolite);
      }

      for(unsigned int i = 0; i < model->getNumReactions(); i++)
      {
          libsbml::Reaction* sbml_reaction = listOfPhysiFBA_reactions->get(i);

          PhysiFBA_reaction* reaction = new PhysiFBA_reaction(sbml_reaction->getId());
          reaction->setName(sbml_reaction->getName());

          libsbml::FbcReactionPlugin* rxnFbc = static_cast<libsbml::FbcReactionPlugin*> (sbml_reaction->getPlugin("fbc"));
          if ( rxnFbc )
          {
              // Getting reaction's upper and lower bounds
              const std::string lbId = rxnFbc->getLowerFluxBound();
              double lb = listOfParameters->get(lbId)->getValue();
              reaction->setLowerBound(lb);

              const std::string ubId = rxnFbc->getUpperFluxBound();
              double ub = listOfParameters->get(ubId)->getValue();
              reaction->setUpperBound(ub);
          }
          int numReactans = sbml_reaction->getNumReactants();
          for(int j = 0; j < numReactans; j++)
          {
              libsbml::SpeciesReference* sbml_species = sbml_reaction->getReactant(j);
              double stoich_coef = -1. * sbml_species->getStoichiometry();

              if ( !this->hasMetabolite(sbml_species->getSpecies()) )
              {
                  PhysiFBA_metabolite* metabolite = new PhysiFBA_metabolite(sbml_species->getSpecies());
                  metabolite->setName(sbml_species->getName());
                  this->addMetabolite(metabolite);
              }
              const PhysiFBA_metabolite* metabolite = this->getMetabolite(sbml_species->getSpecies());
              if (metabolite != nullptr)
                  reaction->addMetabolite(metabolite, stoich_coef);
              else
                  std::cout << "ERROR: PhysiFBA_metabolite " << sbml_species->getSpecies() << " not found" << std::endl;
          }

          int numProducts = sbml_reaction->getNumProducts();
          for(int j = 0; j < numProducts; j++)
          {
              libsbml::SpeciesReference* sbml_species = sbml_reaction->getProduct(j);
              double stoich_coef = sbml_species->getStoichiometry();

              if ( !this->hasMetabolite(sbml_species->getSpecies()) )
              {
                  PhysiFBA_metabolite* metabolite = new PhysiFBA_metabolite(sbml_species->getSpecies());
                  metabolite->setName(sbml_species->getName());
                  addMetabolite(metabolite);
              }
              const PhysiFBA_metabolite* metabolite = this->getMetabolite(sbml_species->getSpecies());
              if (metabolite != nullptr)
                  reaction->addMetabolite(metabolite, stoich_coef);
              else
                  std::cout << "ERROR: PhysiFBA_metabolite " << sbml_species->getSpecies() << " not found" << std::endl;
          }
          this->addReaction(reaction);
      }

      // The following code is intended to extract the objective function from the sbml using
      // the FbcModelPlugin; then the coefficients are assigned to the corresponding reactions
      libsbml::FbcModelPlugin* mplugin = static_cast<libsbml::FbcModelPlugin*>(model->getPlugin("fbc"));
      libsbml::ListOfObjectives* listOfObjectives =  mplugin->getListOfObjectives();
      libsbml::Objective* objective =  mplugin->getObjective(listOfObjectives->getActiveObjective());
      libsbml::ListOfFluxObjectives* listOfFluxObjectives = objective->getListOfFluxObjectives();

      for(unsigned int i=0; i <listOfFluxObjectives->getNumFluxObjectives(); i++)
      {
          libsbml::FluxObjective* fluxObjective = listOfFluxObjectives->get(i);
          std::string rId = fluxObjective->getReaction();
          double objectiveCoefficient = fluxObjective->getCoefficient();
          PhysiFBA_reaction* reaction = this->getReaction(rId);
          reaction->setObjectiveCoefficient(objectiveCoefficient);
      }

      delete document;
  }

  void PhysiFBA_model::initLpModel()
  {

      int n_rows = this->getNumMetabolites();
      int n_cols = this->getNumReactions();

      this->lp_model = new ClpSimplex();

      CoinPackedMatrix matrix;
      matrix.setDimensions(n_rows, 0);

      double* row_lb = new double[n_rows]; //the row lower bounds
      double* row_ub = new double[n_rows]; //the row upper bounds
      double* col_lb = new double[n_cols]; //the column lower bounds
      double* col_ub = new double[n_cols]; //the column upper bounds
      double* objective = new double[n_cols]; //the objective coefficients

      for(int i=0; i< this->getNumMetabolites(); i++)
      {
          row_lb[i] = 0;
          row_ub[i] = 0;
      }
      for(PhysiFBA_reaction* reaction: this->reactions)
      {
          int col_idx = this->reactionsIndexer[reaction->getId()];
          col_lb[col_idx] = reaction->getLowerBound();
          col_ub[col_idx] = reaction->getUpperBound();
          objective[col_idx] = reaction->getObjectiveCoefficient();

          const std::map<const PhysiFBA_metabolite*, double> metabolites = reaction->getMetabolites();
          CoinPackedVector col;
          for(auto it=metabolites.begin(); it!=metabolites.end(); it++)
          {
              const PhysiFBA_metabolite* metabolite = it->first;
              double stoich_coeff = it->second;
              int row_idx = this->metaboliteIndexer[metabolite->getId()];
              col.insert(row_idx, stoich_coeff);
          }
          matrix.appendCol(col);
      }

      this->lp_model->loadProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);
      this->lp_model->setOptimizationDirection(-1);

      delete col_lb;
      delete col_ub;
      delete row_lb;
      delete row_ub;
      delete objective;
  }

  void PhysiFBA_model::writeLp(const char *filename)
  {
    //this->lp_model->writeLp(filename, "lp");
  }

  void PhysiFBA_model::runFBA()
  {
      std::cout << "Before running" << std::endl;
      this->lp_model->primal();
      if ( lp_model->isProvenOptimal() )
      {
          double * columnPrimal = this->lp_model->primalColumnSolution();
          for(PhysiFBA_reaction* reaction: this->reactions)
          {
              int idx = this->reactionsIndexer[reaction->getId()];
              double v = columnPrimal[idx];
              reaction->setFluxValue(v);
          }
      }
      else
      {
          std::cout << "Primal infeasible" << std::endl;
      }
  }

  bool PhysiFBA_model::getSolutionStatus()
  {
      if (this->lp_model)
      {
          return this->lp_model->isProvenOptimal();
      }

      return false;
  }

  float PhysiFBA_model::getObjectiveValue()
  {
      assert(("LP problem must be initialzed first using initLpModel()", this->lp_model != nullptr));
      if (this->lp_model->isProvenOptimal())
          return this->lp_model->getObjValue();
      else
          std::cout << "WARNING: Primal infeasible" << std::endl;
      return 0;
  }


}
