/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "cancer_flux_balance.h"



// declare cell definitions here
// Cell_Definition metabolic_cell;

void create_cell_types( void )
{

	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs

	SeedRandom( parameters.ints( "random_seed" ) );

	// housekeeping

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// turn the default cycle model to live,
	// so it's easier to turn off proliferation

	cell_defaults.phenotype.cycle.sync_to_cycle_model( live );

	cell_defaults.name = "metabolic cell";
	cell_defaults.type = 0;

	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;

	// use default proliferation and death

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0;

	cell_defaults.parameters.o2_proliferation_saturation = 38.0;
	cell_defaults.parameters.o2_reference = 38.0;

	// set oxygen uptake and secretion to zero
	static int oxygen_idx = microenvironment.find_density_index( "oxygen" ); // 0
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_idx] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_idx] = 10;
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_idx] = 38;

	// set the default cell type to no phenotype updates
	cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.volume_update_function = update_cell;

	return;
}

void setup_microenvironment( void )
{

	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "WARNING: overriding from 3-D to 2-D" << std::endl;
		default_microenvironment_options.simulate_2D = true;
	}

	default_microenvironment_options.calculate_gradients = true;

	default_microenvironment_options.track_internalized_substrates_in_each_agent = false;
	default_microenvironment_options.outer_Dirichlet_conditions = true;

	default_microenvironment_options.use_oxygen_as_first_field = true;
	

	initialize_microenvironment();

	return;
}

void setup_tissue( void )
{

	double box_radius = 100.;  // rf. config file
	double xmin = -box_radius;
	double ymin = -box_radius;
	double xmax = box_radius;
	double ymax = box_radius;
	double delta = 4.0;

	Cell* pC;

	pC = create_cell( cell_defaults );
	pC->assign_position( 0., 0.0, 0.0 );
	pC->physifba_model = PhysiFBA::PhysiFBA_default_model;

}

void update_cell(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt ){
  /*
  for(auto itr = PhysiFBA::exchange_flux_density_map.begin();
      itr != PhysiFBA::exchange_flux_density_map.end(); itr++)
    {
      std::string exchange_id = itr->first;
      std::string density_id = itr->second;
      //pCell->physifba_model;
  }
  */

  std::string density_name = "oxygen";

  std::string oxygen_flux_id = PhysiFBA::exchange_flux_density_map[density_name];
  std::string biomass_flux_id = PhysiFBA::exchange_flux_density_map["growth_rate"];

  static int density_idx = microenvironment.find_density_index( density_name );
  double density_value = pCell->nearest_density_vector()[density_idx];
  std::cout << "oxygen around me (" << pCell->ID << ") is: " << density_value << std::endl;
  
  double O2_Km = parameters.doubles("oxygen_Km");
  double O2_Vmax = parameters.doubles("oxygen_Vmax");
  
  double flux_bound = -1 * ( O2_Vmax * density_value) / (density_value + O2_Km);
  pCell->physifba_model.setReactionLowerBound(oxygen_flux_id, flux_bound);

/*
  std::string ex_glc_id = "R_EX_glc__D_e";
  density_name = PhysiFBA::exchange_flux_density_map[ex_glc_id];
  density_idx = microenvironment.find_density_index(density_name);
  density_value = pCell->nearest_density_vector()[density_idx];

  double Glc_Km = 0.0001;
  double Glc_Vmax = 10.;
  double Glc_flux_bound = ( Glc_Vmax * density_value) / (density_value + Glc_Km);
*/


  std::cout << "Running FBA for cell: " << pCell->ID << std::endl;
  pCell->physifba_model.runFBA();
    
  PhysiFBA_reaction* reaction = pCell->physifba_model.getReaction(oxygen_flux_id);
  double oxygen_flux = reaction->getFluxValue();
  std::cout << "Oxygen flux: " << oxygen_flux << std::endl;
  
  PhysiFBA_reaction* biomass_reaction = pCell->physifba_model.getReaction( biomass_flux_id );
  double growht_rate = biomass_reaction->getFluxValue();
  std::cout << "Biomass flux: " << growht_rate << std::endl;

/*
  phenotype.volume.fluid += dt * phenotype.volume.fluid_change_rate *
  	( phenotype.volume.target_fluid_fraction * phenotype.volume.total - phenotype.volume.fluid );

  if( phenotype.volume.fluid < 0.0 )
  { phenotype.volume.fluid = 0.0; }

  phenotype.volume.cytoplasmic_fluid = phenotype.volume.fluid;

  phenotype.volume.cytoplasmic_solid += dt * phenotype.volume.cytoplasmic_biomass_change_rate *
    (phenotype.volume.target_solid_cytoplasmic - phenotype.volume.cytoplasmic_solid );

  if( phenotype.volume.cytoplasmic_solid < 0.0 )
  { phenotype.volume.cytoplasmic_solid = 0.0; }

  phenotype.volume.solid = phenotype.volume.cytoplasmic_solid;

  phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + phenotype.volume.cytoplasmic_fluid;
  phenotype.volume.total = phenotype.volume.cytoplasmic_solid + phenotype.volume.cytoplasmic_fluid;

  // Tell physicell to update the cell radius to the new volume
  phenotype.geometry.update(pCell, phenotype, dt);
*/

}

void setup_default_metabolic_model( void )
{

  // This codes is to use iSIM model
  std::string sbml_fileame = parameters.strings("sbml_model");
  PhysiFBA::PhysiFBA_default_model.readSBMLModel(sbml_fileame.c_str());
  PhysiFBA::PhysiFBA_default_model.initLpModel();

  PhysiFBA::exchange_flux_density_map["glucose"] = "R_E1";
  PhysiFBA::exchange_flux_density_map["lactate"] = "R_E2";
  PhysiFBA::exchange_flux_density_map["oxygen"] = "R_E3";
  PhysiFBA::exchange_flux_density_map["growth_rate"] = "R_R4";

}

void anuclear_volume_model (Cell* pCell, Phenotype& phenotype, double dt)
{

    return;
}

void metabolic_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if cell is dead, don't bother with future phenotype changes.
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	// update the transition rate according to growth rate?
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	//static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" );
	//phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
	return;
}



std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring

	std::vector<std::string> output = false_cell_coloring_cytometry(pCell);
	output[0] = "red";
	output[1] = "red";
	output[2] = "red";

	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		 output[0] = "black";
		 output[2] = "black";
	}

	return output;
}
