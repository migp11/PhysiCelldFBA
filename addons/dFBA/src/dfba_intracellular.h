#ifndef _dFBA_Intracellular_h_
#define _dFBA_Intracellular_h_

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>   // for setw

#include <coin/CoinPackedMatrix.hpp>
#include <coin/CoinPackedVector.hpp>
#include <coin/ClpSimplex.hpp>

#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"

#include "dfba_Model.h"

using namespace std;

namespace PhysiCelldFBA {

static std::string PhysiCelldFBA_Version = "0.0.1"; 

static float hours_to_minutes = 1/60;
static float PI = PhysiCell::PhysiCell_constants::pi;

struct KineticParam
{
	string name;
	string untis;
	float value;
};

struct ExchangeFluxData
{
	string density_name;
	string fba_flux_id;
	int density_index;
	KineticParam Km;
	KineticParam Vmax;
};

void update_dfba_inputs( PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt );

void update_dfba_outputs( PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt );



class dFBAIntracellular : public PhysiCell::Intracellular 
{
 private:
 public:

 	std::string sbml_filename;
	double reference_volume = 1.0;
	double cell_density = 1.04;
	double max_growth_rate = 0;
	double current_growth_rate = 0;
	double next_dfba_run = 0;
	dFBAModel model;

	ClpSimplex problem;
	CoinMessageHandler* handler;
	bool is_initialized = false;

	/** \brief map between metabolites' ids and metabolites' references **/
	std::map<std::string, int> metaboliteIndexer;

	/** \brief map between reaction IDs and reaction references */
	std::map< std::string, int> reactionsIndexer;

	/** \brief map between density IDs and exchange reactions */
	std::map<std::string, ExchangeFluxData> substrate_exchanges;

    dFBAIntracellular();

	dFBAIntracellular(pugi::xml_node& node);
	
	dFBAIntracellular(dFBAIntracellular* copy);

		// rwh: review this
	Intracellular* clone(){
		return static_cast<Intracellular*>(new dFBAIntracellular(this));
	}

	Intracellular* getIntracellularModel() 
	{
		return static_cast<Intracellular*>(this);
	}
	
    // ================  generic  ================
	// This function parse the xml cell definition
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
	// This function checks if it's time to update the model
	bool need_update() { return PhysiCell::PhysiCell_globals.current_time >= this->next_dfba_run; }

	// This function deals with inheritance from mother to daughter cells
	void inherit(PhysiCell::Cell* cell){

	};

	// Get value for model parameter
	double get_parameter_value(std::string name) { return 0;}
	
	// Set value for model parameter
	void set_parameter_value(std::string name, double value) { return; }

	std::string get_state(){ return ""; }
	
	void display(std::ostream& os){ return; }
	
//	~Intracellular();
	

	
	int parse_transport_model(pugi::xml_node& node);
	void parse_growth_model(pugi::xml_node& node);
	void initLpSolver();

		
	void update();
	void update(PhysiCell::Cell* cell, PhysiCell::Phenotype& phenotype, double dt){
		return;
	};
	
	
	// libroadrunner specifics
		
	// for now, define dummy methods for these in the abstract parent class
	
	// This function initialize the model, needs to be called on each cell once created
	void start();
	
	
	// static void save_PhysiBoSS(std::string path, std::string index);
	static void save_dFBA(std::string path, std::string index);

	void update_volume(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double growth_rate, double dt);
	void standard_update_cell_volume(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double growth_rate, double dt);



	// unneeded for this type

    // ================  specific to "maboss" ================
	bool has_variable(std::string name) { return false; }
	bool get_boolean_variable_value(std::string name) { return false; }
	void set_boolean_variable_value(std::string name, bool value) {	}
	void print_current_nodes(){	}
	

    // ================  specific to "roadrunner" ================
    int update_phenotype_parameters(PhysiCell::Phenotype& phenotype) {return 0; }
    int validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype) {return 0; }
    int validate_SBML_species() {return 0; }
    int create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype) {return 0; }

};

}

#endif
