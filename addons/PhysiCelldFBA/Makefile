# C++ Compiler command
CXX = g++

COINUTILS_HAS_GLPK = false

BUILD_DIR = ./build
SRC_DIR = ./src
BIN_DIR = ./bin
CUR_DIR = $(shell pwd)

# Folders to look in
vpath %.o     $(BUILD_DIR)
vpath %.cpp   $(SRC_DIR) 
vpath %.h     $(SRC_DIR) 


# Path to the conda env where coincbc and libsbml-plus-packaes were installed
# conda install -c conda-forge/label/gcc7 coincbc
# conda install -c vincent-noel libsbml-plus-package

DEPS_DIR = /home/mponce/.conda/envs/pc-dfba

CPPFLAGS  = -I$(DEPS_DIR)/include
LDFLAGS   = -L$(DEPS_DIR)/lib

# CLP_DIR = $(CUR_DIR)/ext/coin-or
# LIBSBML_DIR = $(CUR_DIR)/ext/libsbml
# Icnclude and folder libs
# CPPFLAGS  += -I$(CLP_DIR)/include
# CPPFLAGS  += -I$(LIBSBML_DIR)/include
# LDFLAGS   += -L$(CLP_DIR)/lib
# LDFLAGS   += -L$(LIBSBML_DIR)/lib64

# Libraries (other: -lcoinasl -lcoinglpk -lcoinmumps)
LIBS      += -lsbml
LIBS      += -lClp -lCoinUtils 
LIBS      += -lstdc++ -ldl -lm

# Additionals flags for debugging -g3 -fno-rtti -ggdb

ARCH := native
CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -pipe 

COMPILE_COMMAND := $(CXX) $(CFLAGS) 

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o \
PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_basic_signaling.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o PhysiCell_geometry.o

pugixml_OBJECTS := pugixml.o

PC_OBJECTS := $(BioFVM_OBJECTS) $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)

FBA_OBJECTS := FBA_model.o FBA_reaction.o FBA_metabolite.o
INTRACELLULAR := dfba_intracellular.o
ALL_OBJECTS := $(FBA_OBJECTS) $(INTRACELLULAR)

make_lib: $(ALL_OBJECTS) $(PC_OBJECTS)
programs: testSBML testClp testFBA testIntracellular
all: make_lib programs

build:
	@mkdir -p $@

$(ALL_OBJECTS): %.o: %.cpp %.h | build 
	$(COMPILE_COMMAND) $(CPPFLAGS) $(LDFLAGS) -c $< -o $(BUILD_DIR)/$(@F) $(LIBS)

# create test just for libSBML
testSBML: ./test/testSBML.cpp
	$(COMPILE_COMMAND) -I./usr/include -L./usr/lib -o $(BIN_DIR)/testSBML ./test/testSBML.cpp -lsbml-static
# $(COMPILE_COMMAND) $(CPPFLAGS) $(LDFLAGS) -o $(BIN_DIR)/testSBML ./test/testSBML.cpp -lsbml

# create test just for Clp
testClp: ./test/testClp.cpp
	$(COMPILE_COMMAND) $(CPPFLAGS) $(LDFLAGS) -o $(BIN_DIR)/testClp ./test/testClp.cpp -lClp -lCoinUtils

testFBA: ./test/testFBA.cpp $(FBA_OBJECTS)
	$(COMPILE_COMMAND) $(CPPFLAGS) $(LDFLAGS) -o $(BIN_DIR)/testFBA ./test/testFBA.cpp $(BUILD_DIR)/FBA_*.o $(LIBS)

testIntracellular: ./test/testIntracellular.cpp $(ALL_OBJECTS) $(PC_OBJECTS)
	$(COMPILE_COMMAND) $(CPPFLAGS) $(LDFLAGS) -o $(BIN_DIR)/testIntracellular ./test/testIntracellular.cpp $(BUILD_DIR)/*.o $(LIBS)

clean:
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BIN_DIR)/* 



MODULES := ../../modules
CORE := ../../core
BIOFVM := ../../BioFVM


PhysiCell_phenotype.o: $(CORE)/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_phenotype.cpp -o $(BUILD_DIR)/$@
	
PhysiCell_digital_cell_line.o: $(CORE)/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_digital_cell_line.cpp -o $(BUILD_DIR)/$@

PhysiCell_cell.o: $(CORE)/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_cell.cpp  -o $(BUILD_DIR)/$@

PhysiCell_cell_container.o: $(CORE)/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_cell_container.cpp  -o $(BUILD_DIR)/$@
	
PhysiCell_standard_models.o: $(CORE)/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_standard_models.cpp  -o $(BUILD_DIR)/$@
	
PhysiCell_utilities.o: $(CORE)/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_utilities.cpp  -o $(BUILD_DIR)/$@
	
PhysiCell_custom.o: $(CORE)/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_custom.cpp -o $(BUILD_DIR)/$@
	
PhysiCell_constants.o: $(CORE)/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_constants.cpp  -o $(BUILD_DIR)/$@
	
# BioFVM core components (needed by PhysiCell)
	
BioFVM_vector.o: $(BIOFVM)/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_vector.cpp -o $(BUILD_DIR)/$@

BioFVM_agent_container.o: $(BIOFVM)/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_agent_container.cpp -o $(BUILD_DIR)/$@
	
BioFVM_mesh.o: $(BIOFVM)/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_mesh.cpp 

BioFVM_microenvironment.o: $(BIOFVM)/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_microenvironment.cpp -o $(BUILD_DIR)/$@

BioFVM_solvers.o: $(BIOFVM)/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_solvers.cpp -o $(BUILD_DIR)/$@

BioFVM_utilities.o: $(BIOFVM)/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_utilities.cpp -o $(BUILD_DIR)/$@
	
BioFVM_basic_agent.o: $(BIOFVM)/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_basic_agent.cpp -o $(BUILD_DIR)/$@
	
BioFVM_matlab.o: $(BIOFVM)/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_matlab.cpp -o $(BUILD_DIR)/$@

BioFVM_MultiCellDS.o: $(BIOFVM)/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/BioFVM_MultiCellDS.cpp -o $(BUILD_DIR)/$@
	
pugixml.o: $(BIOFVM)/pugixml.cpp
	$(COMPILE_COMMAND) -c $(BIOFVM)/pugixml.cpp -o $(BUILD_DIR)/$@
	
# standard PhysiCell modules

PhysiCell_SVG.o: $(MODULES)/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c $(MODULES)/PhysiCell_SVG.cpp -o $(BUILD_DIR)/$@
    
PhysiCell_pathology.o: $(MODULES)/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c $(MODULES)/PhysiCell_pathology.cpp -o $(BUILD_DIR)/$@

PhysiCell_MultiCellDS.o: $(MODULES)/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c $(MODULES)/PhysiCell_MultiCellDS.cpp -o $(BUILD_DIR)/$@

PhysiCell_various_outputs.o: $(MODULES)/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c $(MODULES)/PhysiCell_various_outputs.cpp -o $(BUILD_DIR)/$@
	
PhysiCell_pugixml.o: $(MODULES)/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c $(MODULES)/PhysiCell_pugixml.cpp -o $(BUILD_DIR)/$@
	
PhysiCell_settings.o: $(MODULES)/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c $(MODULES)/PhysiCell_settings.cpp	-o $(BUILD_DIR)/$@
	
PhysiCell_basic_signaling.o: $(CORE)/PhysiCell_basic_signaling.cpp
	$(COMPILE_COMMAND) -c $(CORE)/PhysiCell_basic_signaling.cpp -o $(BUILD_DIR)/$@

PhysiCell_geometry.o: $(MODULES)/PhysiCell_geometry.cpp
	$(COMPILE_COMMAND) -c $(MODULES)/PhysiCell_geometry.cpp -o $(BUILD_DIR)/$@
