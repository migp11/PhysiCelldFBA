# C++ Compiler command
CXX = g++

COINUTILS_HAS_GLPK = false

BUILD_DIR = ./build
SRC_DIR = ./src
BIN_DIR = ./bin
CUR_DIR = $(shell pwd)
CLP_DIR = $(CUR_DIR)/ext/coin-or
LIBSBML_DIR = $(CUR_DIR)/ext/libsbml

# Folders to look in
vpath %.o     $(BUILD_DIR)
vpath %.cpp   $(SRC_DIR) 
vpath %.h     $(SRC_DIR) 

# Libraries (other: -lcoinasl -lcoinglpk -lcoinmumps)
LIBS      += -lsbml 
LIBS      += -lClp -lCoinUtils 
LIBS      += -lstdc++ -ldl -lm

# Icnclude and folder libs
CPPFLAGS  += -I$(CLP_DIR)/include
CPPFLAGS  += -I$(LIBSBML_DIR)/include
LDFLAGS   += -L$(CLP_DIR)/lib
LDFLAGS   += -L$(LIBSBML_DIR)/lib64

# Additionals flags for debugging -g3 -fno-rtti -ggdb
CXXFLAGS= -Wall -std=c++11 -pipe 

COMPILE_COMMAND := $(CXX) $(CXXFLAGS) 

FBA_OBJECTS := FBA_model.o FBA_reaction.o FBA_metabolite.o
ALL_OBJECTS := $(FBA_OBJECTS)

make_lib: $(ALL_OBJECTS)
programs: testSBML testClp testFBA
all: make_lib programs

build:
	@mkdir -p $@

$(ALL_OBJECTS): %.o: %.cpp %.h | build 
	$(COMPILE_COMMAND) $(CPPFLAGS) $(LDFLAGS) -c $< -o $(BUILD_DIR)/$(@F) $(LIBS)

# create test just for libSBML
testSBML: ./test/testSBML.cpp
	$(COMPILE_COMMAND) -o $(BIN_DIR)/testSBML ./test/testSBML.cpp -lsbml
# $(COMPILE_COMMAND) -L$(LIBSBML_DIR)/usr/lib64 -I$(LIBSBML_DIR)/usr/include -o $(BIN_DIR)/testSBML ./test/testSBML.cpp -lsbml

# create test just for Clp
testClp: ./test/testClp.cpp
	$(COMPILE_COMMAND) -L$(CLP_DIR)/lib -I$(CLP_DIR)/include -o $(BIN_DIR)/testClp ./test/testClp.cpp -lClp -lCoinUtils

testFBA: ./test/testFBA.cpp $(FBA_OBJECTS)
	$(COMPILE_COMMAND) $(CPPFLAGS) $(LDFLAGS) -o $(BIN_DIR)/testFBA ./test/testFBA.cpp $(BUILD_DIR)/*.o $(LIBS)

clean:
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BIN_DIR)/* 
