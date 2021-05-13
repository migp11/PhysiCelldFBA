# C++ Compiler command
CXX = g++

COINUTILS_HAS_GLPK = false

BUILD_DIR = build
SRC_DIR = src
BIN_DIR = bin
CUR_DIR = $(shell pwd)

CLP_DIR = $(CUR_DIR)/coin-or

INC := -I/usr/local/include -I$(CLP_DIR)/include
# LIBS := -L/usr/local/lib -L$(CLP_DIR)/lib -lClp -lCoinUtils -lsbml -ldl -lm -lcoinasl -lcoinglpk -lcoinmumps
LIBS := -L/usr/local/lib -L$(CLP_DIR)/lib -lClp -lCoinUtils -lsbml -ldl -lm -lcoinasl -lcoinglpk -lcoinmumps

# CFLAGS will be the options passed to the compiler.
# CFLAGS= -c -Wall -lstdc++fs -std=c++17 -g3 -fno-rtti -ggdb

#---------- Folders to look in
vpath %.o     $(BUILD_DIR)
vpath %.cpp   $(SRC_DIR) 
vpath %.h     $(SRC_DIR) 

# Additionals flags for debugging -g3 -fno-rtti -ggdb
# -O0 -ggdb
CXXFLAGS= -O3 -Wall -std=c++11 -pipe 

COMPILE_COMMAND := $(CXX) $(CXXFLAGS) 

FBA_OBJECTS := FBA_model.o FBA_reaction.o FBA_metabolite.o
ALL_OBJECTS := $(FBA_OBJECTS)
MAIN1 := ./test/main_test_fba.cpp
MAIN2 := ./test/main_test_sbml.cpp


make_lib: $(ALL_OBJECTS)
all: make_lib test_fba test_sbml


build:
	@mkdir -p $@

$(ALL_OBJECTS): %.o: %.cpp %.h | build 
	$(COMPILE_COMMAND) $(CXXFLAGS) $(INC) $(LIBS) -c $< -o $(BUILD_DIR)/$(@F)

test_fba: $(MAIN1)
	$(COMPILE_COMMAND) $(INC) -o $(BIN_DIR)/test_fba ./test/main_test_fba.cpp $(BUILD_DIR)/*.o $(LIBS)

test_sbml: $(BIN_DIR) $(MAIN2) 
	$(COMPILE_COMMAND) $(INC) -o $(BIN_DIR)/test_sbml ./test/main_test_sbml.cpp $(BUILD_DIR)/*.o $(LIBS)

clean:
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BIN_DIR)/* 
