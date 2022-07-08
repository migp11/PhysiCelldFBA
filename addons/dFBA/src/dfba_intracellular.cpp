#include "dfba_intracellular.h"
#include <sstream>
#include <iostream>

dFBAIntracellular::dFBAIntracellular() : Intracellular()
{
    this->intracellular_type = "dfba";
}

dFBAIntracellular::dFBAIntracellular(pugi::xml_node& node)
{
    intracellular_type = "dfba";
	this->initialize_intracellular_from_pugixml(node);
}

dFBAIntracellular::dFBAIntracellular(dFBAIntracellular* copy) 
{
    intracellular_type = copy->intracellular_type;
    sbml_filename = copy->sbml_filename;
    parameters = copy->parameters;
    fba_model.initFBAmodel(copy->sbml_filename.c_str());
}


void dFBAIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{

    pugi::xml_node node_sbml = node.child( "sbml_filename" );

	if ( node_sbml )
	{ 
        this->sbml_filename = PhysiCell::xml_get_my_string_value (node_sbml);
        
    }
    else
    {
        std::cout << "Error: attempted get sbml_filename but not foun. " << std::endl;
        std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }
	
    pugi::xml_node node_exchange = node.child( "exchange" );
    
	while( node_exchange )
	{
        exchange_data ex_struct;
        kinetic_parm Km;
        kinetic_parm Vmax;
        

		string density_name = node_exchange.attribute( "substrate" ).value(); 
        int density_index = microenvironment.find_density_index( density_name ); 
        std::string actual_name = microenvironment.density_names[ density_index ]; 
			
        // error check 
        if( std::strcmp( density_name.c_str() , actual_name.c_str() ) != 0 )
        {
            std::cout << "Error: attempted to set secretion/uptake/export for \""; 
            std::cout << density_name << "\", which was not found in the microenvironment." << std::endl;
            std::cout << "Please double-check your substrate name in the config file." << std::endl;
            std::cout << std::endl;
            exit(-1); 
        }
        
        pugi::xml_node node_fba_flux = node_exchange.child( "fba_flux" ); 
		if( node_fba_flux )
		{  
            ex_struct.fba_flux_id = PhysiCell::xml_get_my_string_value(node_fba_flux);
        }
        else {
            std::cout << "Error: attempted get fba_flux node for "; 
            std::cout << ex_struct.density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }

        pugi::xml_node node_Km = node_exchange.child( "Km" ); 
		if( node_Km )
		{
            Km.name = "Km";
            Km.untis = node_Km.attribute("units").value();
            Km.value = PhysiCell::xml_get_my_double_value(node_Km);
            
        }
        else {
            std::cout << "Error: attempted get Km node for "; 
            std::cout << density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }

        pugi::xml_node node_Vmax = node_exchange.child( "Vmax" ); 
		if( node_Vmax )
		{
            Vmax.name = "Vmax";
            Vmax.untis = node_Vmax.attribute("units").value();
            Vmax.value = PhysiCell::xml_get_my_double_value(node_Vmax);
            
        }
        else {
            std::cout << "Error: attempted get Vmax node for "; 
            std::cout << density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }
		
        ex_struct.density_name = density_name;
        ex_struct.density_index = density_index;
        ex_struct.Km = Km;
        ex_struct.Vmax = Vmax;

        this->substrate_exchanges[density_name] = ex_struct;
		node_exchange = node_exchange.next_sibling( "exchange" ); 
	}

    std::cout << "Loaing SBML model from: " << this->sbml_filename << std::endl;
    this->fba_model.readSBMLModel(this->sbml_filename.c_str());
    this->fba_model.initLpModel();
    this->fba_model.runFBA();

}

void dFBAIntracellular::dFBAIntracellular::start()
{
    return;
}

bool dFBAIntracellular::dFBAIntracellular::need_update()
{
    return PhysiCell::PhysiCell_globals.current_time >= this->next_model_run;
}

void dFBAIntracellular::update(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt)
{
    /*
        Steps for the update
        1- Update exchange fluxes lower bound using current concentrations
        2- run FBA and check growth threshold
        3- update the cell volumne using the growth rate from FBA
        4- rescale exchange fluxes from the dfba model and use them to update the net_export_rates
        5- remove the internalized substrates if needed

        mM /gDW cell / h <-- fluxes uints
        Biomass: 

    */
    
    static float pi = PhysiCell::PhysiCell_constants::pi;
    bool debug = false;

    // HeLa cell mass 2.3 ng

    // cell volume fL (um³)
    // mM: mmol / L = 10⁻³ mol / 1e¹⁵ um³ = 10⁻¹⁵mol / um³ = amol /  um³
    // cell_density --> 1.04 g/ml = 1.04 ug/nL = 1.04 ng / pL = 1.04 pg / fL = pg / um³
    // cell.mass.total (ng) = cell.volume (um³) * cell.density (pg / um³) = pg
    // ~2500 (um³) * ~ 1.04 pg / (um³) =~ 2600 (pg)
    // gDW cell (cell.volume.total * cell_density) = mass.total
    // cell.mass.solid = cell.mass.total * (1-fluid_frac)

    static float cell_density = 1.04;                             // pg / um³
    float solid_fraction = 1 - phenotype.volume.fluid_fraction;   // unitless ~30%
    float solid_volume = phenotype.volume.total * solid_fraction; // um³
    float cell_dry_weight = solid_volume * cell_density ;         // um³ * pg / um³ = pg

    // r = (3V / 4π))^1/3
    float radius = cbrt( (3./4. * pi * phenotype.volume.total) ); // um
    float cell_surface = 4 * pi * pow(radius, 2);                 // um

    std::vector<double> density_vector = pCell->nearest_density_vector(); // 

    // Setp 1 - Update exchange fluxes lower
    map<std::string, exchange_data>::iterator it;
    for(it = this->substrate_exchanges.begin(); it != this->substrate_exchanges.end(); it++)
    {
        std::string substrate_name = it->first;
        exchange_data ex_strut = it->second;

        double Vmax = ex_strut.Vmax.value;
        double Km = ex_strut.Km.value;
        
        // geting the amount of substrate
        double substrate_conc = density_vector[ex_strut.density_index];;

        // Standard FBA fluxes units:  mmol /  gDW cell / h
        // Km: mM = mmol/L
        // sdt: surf_density_of_transport
        // kcat: catalytic constant
        // Vmax = kcat * sdt * cell_surface
        
        // using irreversible Michaelis Menten kinetics to estimate the flux bound
        double flux_bound = cell_dry_weight * (Vmax * substrate_conc) / (Km + substrate_conc); // should be calculated from density
        // Change sign to use as lower bound of the exchange flux
        flux_bound *= -1;
        // Updateing the lower bound of the corresponding exchange flux
        this->fba_model.setReactionLowerBound(ex_strut.fba_flux_id, flux_bound);    
        
        if ( debug ) {
            std::cout << " - [" << substrate_name << "] = " << substrate_conc;
            std::cout << " ==> " << ex_strut.fba_flux_id << " = " << flux_bound << std::endl;
        }
    }

    static float hours_to_minutes = 1/60;

    // STEP 2 - run FBA
    this->fba_model.runFBA();

    // STEP 3 - Update cell volumen using growth rate (first rscale growth rate to 1/min)
    float growth_rate = this->fba_model.getObjectiveValue();
    growth_rate = growth_rate * hours_to_minutes; // growth_rate 1/h -> 1/min
    
    // V(t+1) = V(t) + V(t) * mu * dt = V(t) * (1 + u * dt) 
    float volume_increase_ratio = 1 + growth_rate * dt;
    phenotype.volume.multiply_by_ratio( volume_increase_ratio );
    phenotype.geometry.update(pCell, phenotype, dt);

    
    // STEPS 4-5 - Update net_export_rates for the different densities
    for(it = this->substrate_exchanges.begin(); it != this->substrate_exchanges.end(); it++)
    {
        // Retrive the exchange flux and its corresponding density
        std::string substrate_name = it->first;
        exchange_data ex_strut = it->second;
        
        int density_index = ex_strut.density_index;
        std::string fba_flux_id = ex_strut.fba_flux_id;
        
        FBA_reaction* exchange_flux = this->fba_model.getReaction(fba_flux_id);
        double flux_value = exchange_flux->getFluxValue(); // mmol/gDW/h

        // Rescaling FBA exchanges flux into net_export_rates
        // Net export rates are expressed in substance/time
        // mmol/gDW/h --> mmol/min
        // flux_value / 60 * cell_dry_weight  = pico*mmol/min = fmol/min
        // fmol/min * 10e3 =  pmol/min
        double net_export_rate = flux_value * hours_to_minutes * cell_dry_weight * 1e3;

        phenotype.secretion.net_export_rates[density_index] = net_export_rate;
        
        if (default_microenvironment_options.track_internalized_substrates_in_each_agent)
        {
            phenotype.molecular.internalized_total_substrates[density_index] = 0;
        }
    }
}

int dFBAIntracellular::update_phenotype_parameters(PhysiCell::Phenotype& phenotype)
{

    return 0;
}




void dFBAIntracellular::save_dFBA(std::string path, std::string index) 
{
	
}