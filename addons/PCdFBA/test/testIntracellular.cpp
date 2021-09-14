#include <iostream>
#include <vector>

#include "../src/dfba_intracellular.h"

using namespace std;


int main (int argc, const char *argv[])
{

    const char *filename;
	bool XML_status = false;

    pugi::xml_document physicell_config_doc; 	
    pugi::xml_node physicell_config_root; 
	
	if( argc > 1 )
	{
        filename = argv[1];
		XML_status = PhysiCell::load_PhysiCell_config_file( filename ); 
		
	}
	else
	{
		cout << "Error no physicell xml config file provided";
        exit(-1);
	}
	if( !XML_status )
	{ exit(-1); }

    cout << "Using config file " << filename << " ... " << std::endl ; 

	
    
    pugi::xml_parse_result result = physicell_config_doc.load_file( filename  );
	
	if( result.status != pugi::xml_parse_status::status_ok )
	{
		std::cout << "Error loading " << filename << "!" << std::endl; 
		exit(-1);
	}
	
	physicell_config_root = physicell_config_doc.child("PhysiCell_settings");
    
	pugi::xml_node node = PhysiCell::xml_find_node(physicell_config_root, "cell_definitions");
	cout << node << endl;
	node = PhysiCell::xml_find_node(node, "cell_definition");
	cout << node << endl;
	node = PhysiCell::xml_find_node(node, "phenotype");
	cout << node << endl;
	node = PhysiCell::xml_find_node(node, "intracellular");
    cout << node << endl;

	pugi::xml_node node_sbml = node.child( "sbml_filename" );
	string sbml_filename = PhysiCell::xml_get_my_string_value (node_sbml);

	cout << sbml_filename << endl;
	dFBAIntracellular* model = new dFBAIntracellular(node);
	cout << "ok" << endl;
	delete model;
	return 0;
}