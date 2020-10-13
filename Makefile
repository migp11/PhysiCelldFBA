VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := ecoli-dfba

CC := g++
# CC := g++-mp-7 # typical macports compiler name
# CC := g++-7 # typical homebrew compiler name

# Check for environment definitions of compiler
# e.g., on CC = g++-7 on OSX
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

ARCH := native # best auto-tuning
# ARCH := core2 # a reasonably safe default for most CPUs since 2007
# ARCH := corei7
# ARCH := corei7-avx # earlier i7
# ARCH := core-avx-i # i7 ivy bridge or newer
# ARCH := core-avx2 # i7 with Haswell or newer
# ARCH := nehalem
# ARCH := westmere
# ARCH := sandybridge # circa 2011
# ARCH := ivybridge   # circa 2012
# ARCH := haswell     # circa 2013
# ARCH := broadwell   # circa 2014
# ARCH := skylake     # circa 2015
# ARCH := bonnell
# ARCH := silvermont
# ARCH := skylake-avx512
# ARCH := nocona #64-bit pentium 4 or later


# For debugging with gdb just add: -g3 -ggdb -O3
# -DCLP_BUILD is an additional flags to link againt clp
# CFLAGS := -march=$(ARCH) -ggdb -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -DCLP_BUILD
CFLAGS := -march=$(ARCH) -ggdb -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -DCLP_BUILD


# Libraries need: SBML and coin-Clp
INC := -I/usr/local/include
LIBS := -L/usr/local/lib -ldl -lm -lClpSolver -lClp -lCoinUtils -lsbml

COMPILE_COMMAND := $(CC) $(CFLAGS)

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o

# put your custom objects here (they should be in the custom_modules directory)
FBA_addon_OBJECTS := FBA_metabolite.o FBA_reaction.o FBA_model.o
PhysiCell_custom_module_OBJECTS := ecoli_acetic_switch.o

pugixml_OBJECTS := pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS)  $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)
<<<<<<< HEAD
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS)

EXAMPLES := ./examples/PhysiCell_test_mechanics_1.cpp ./examples/PhysiCell_test_mechanics_2.cpp \
 ./examples/PhysiCell_test_DCIS.cpp ./examples/PhysiCell_test_HDS.cpp \
 ./examples/PhysiCell_test_cell_cycle.cpp ./examples/PhysiCell_test_volume.cpp 

all: 
	make heterogeneity-sample
	make 

# sample projects 	
list-projects:
	@echo "Sample projects: template2D template3D biorobots-sample cancer-biorobots-sample heterogeneity-sample"
	@echo "                 cancer-immune-sample virus-macrophage-sample ecoli-acetic-switch"

ecoli-acetic-switch:
	cp ./sample_projects/ecoli_acetic_switch/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/ecoli_acetic_switch/main_ecoli_acetic_switch.cpp ./main.cpp
	cp Makefile Makefile-backup
	cp ./sample_projects/ecoli_acetic_switch/Makefile ./
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml
	cp ./sample_projects/ecoli_acetic_switch/config/* ./config/

template2D: 
	cp ./sample_projects/template2D/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/template2D/main-2D.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/template2D/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/template2D/config/* ./config/
	
template3D: 	
	cp ./sample_projects/template3D/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/template3D/main-3D.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/template3D/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/template3D/config/* ./config/
	
# sample projects 

biorobots-sample:
	cp ./sample_projects/biorobots/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/biorobots/main-biorobots.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/biorobots/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/biorobots/config/* ./config/
	
cancer-biorobots-sample:
	cp ./sample_projects/cancer_biorobots/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/cancer_biorobots/main-cancer_biorobots.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/cancer_biorobots/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/cancer_biorobots/config/* ./config/
	
heterogeneity-sample:
	cp ./sample_projects/heterogeneity/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/heterogeneity/main-heterogeneity.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/heterogeneity/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/heterogeneity/config/* ./config/
	
cancer-immune-sample:
	cp ./sample_projects/cancer_immune/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/cancer_immune/main-cancer_immune_3D.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/cancer_immune/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/cancer_immune/config/* ./config/
	
virus-macrophage-sample:
	cp ./sample_projects/virus_macrophage/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/virus_macrophage/main-virus_macrophage.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/virus_macrophage/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/virus_macrophage/config/* ./config/
	
beta-testing:
	cp ./sample_projects/beta_testing/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/beta_testing/main-beta.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/beta_testing/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/beta_testing/config/* ./config/
	
# early examples for convergence testing 

physicell_test_mech1: $(PhysiCell_OBJECTS) ./examples/PhysiCell_test_mechanics_1.cpp 
	$(COMPILE_COMMAND) -o test_mech1 $(PhysiCell_OBJECTS) ./examples/PhysiCell_test_mechanics_1.cpp
=======
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS) $(FBA_addon_OBJECTS) 
>>>>>>> 5fb0001207b22eb032ece6ef408df8fb9f1df5b4

# compile the project

all: main.cpp $(ALL_OBJECTS)
	$(COMPILE_COMMAND) $(INC) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp $(LIBS)

# PhysiCell core components

PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp

PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ./core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell.cpp

PhysiCell_cell_container.o: ./core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_container.cpp

PhysiCell_standard_models.o: ./core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_standard_models.cpp

PhysiCell_utilities.o: ./core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_utilities.cpp

PhysiCell_custom.o: ./core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_custom.cpp
	
PhysiCell_constants.o: ./core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_constants.cpp 

# BioFVM core components (needed by PhysiCell)

BioFVM_vector.o: ./BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_vector.cpp

BioFVM_agent_container.o: ./BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_agent_container.cpp

BioFVM_mesh.o: ./BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_mesh.cpp

BioFVM_microenvironment.o: ./BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_microenvironment.cpp

BioFVM_solvers.o: ./BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_solvers.cpp

BioFVM_utilities.o: ./BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_utilities.cpp

BioFVM_basic_agent.o: ./BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_basic_agent.cpp

BioFVM_matlab.o: ./BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ./BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_MultiCellDS.cpp

pugixml.o: ./BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/pugixml.cpp

# standard PhysiCell modules

PhysiCell_SVG.o: ./modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ./modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ./modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ./modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_various_outputs.cpp

PhysiCell_pugixml.o: ./modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pugixml.cpp

PhysiCell_settings.o: ./modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_settings.cpp

# user-defined PhysiCell modules

# PhysiCelldFBA addon modules

FBA_model.o: ./addons/PhysiCelldFBA/FBA_model.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ./addons/PhysiCelldFBA/FBA_model.cpp -o ./FBA_model.o

FBA_reaction.o: ./addons/PhysiCelldFBA/FBA_reaction.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ./addons/PhysiCelldFBA/FBA_reaction.cpp -o ./FBA_reaction.o

FBA_metabolite.o: ./addons/PhysiCelldFBA/FBA_metabolite.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ./addons/PhysiCelldFBA/FBA_metabolite.cpp -o ./FBA_metabolite.o


ecoli_acetic_switch.o: ./custom_modules/ecoli_acetic_switch.cpp
	$(COMPILE_COMMAND) -c ./custom_modules/ecoli_acetic_switch.cpp


# cleanup

reset:
	rm -f *.cpp
	cp ./sample_projects/Makefile-default Makefile
	rm -f ./custom_modules/*
	touch ./custom_modules/empty.txt
	touch ALL_CITATIONS.txt
	rm ALL_CITATIONS.txt
	cp ./config/PhysiCell_settings-backup.xml ./config/PhysiCell_settings.xml
	touch ./config/empty.txt

clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME)*

data-cleanup:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
<<<<<<< HEAD
	rm -rf ./output
	mkdir ./output
	touch ./output/empty.txt
	
# archival 
=======
	rm -f ./output/*
	touch ./config/empty.txt

# archival
>>>>>>> 5fb0001207b22eb032ece6ef408df8fb9f1df5b4

checkpoint: 
	zip -r $$(date +%b_%d_%Y_%H%M).zip Makefile *.cpp *.h config/*.xml custom_modules/* 

zip:
	zip -r latest.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/*
	cp latest.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest.zip VERSION_$(VERSION).zip
	mv *.zip archives/

tar:
	tar --ignore-failed-read -czf latest.tar Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/*
	cp latest.tar $$(date +%b_%d_%Y_%H%M).tar
	cp latest.tar VERSION_$(VERSION).tar
	mv *.tar archives/

unzip:
	cp ./archives/latest.zip .
	unzip latest.zip

untar:
	cp ./archives/latest.tar .
	tar -xzf latest.tar
