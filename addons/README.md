# PhysiCellFBA
A PhysiCell extenstion to integrate Flux Balance models into cell agents

I've included the sources of the both packages used:
* coin-Clp
* libsbml

So the first step is building both libraries

## On Linux

### Building and installing the coin-Clp librabry:

Some dependencies include Fortran compiler and the libraries BLAS and LAPACK

```
sudo apt install gfortran

sudo apt install libblas-dev liblapack-dev

```

The easiet way to build coin-Clp is by using brew.
So the first step is downloading the coinbrew scrip grom github:
```
git clone https://github.com/coin-or/coinbrew.git
cd coinbrew/
```
The next steps are fetching, building the libriary and its dependencies
```
./coinbrew fetch --main-proj=Clp
mkdir build
./coinbrew build --main-proj=Clp --build-dir=./build --prefix=/usr/local --no-prompt
```
Finally, to install Clp library just run the follow command line as sudo
```
sudo ./coinbrew install --main-proj=Clp
```
That's it! :+1:


### Building and installing the libSBML librabry:

Some dependencies:
```
sudo apt install libxml2-dev 
sudo apt install libbz2-dev
```

To get the current latest version (5.18) source code core plus packages follow the link
(https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/)
The file name will look somthing like: libSBML-5.18.0-core-plus-packages-src.tar.gz

Then just follow the instructions (I recommend using cmake-gui)
**Important** options required:
* ENABLE_FBC
* ENABLE_L3V2EXTENDEDMATH
* WITH_CPP_NAMESPACE


If tha both libraries were correctly installed at /usr/local
you will be able to compile and link the runFBA code

To compile the two binary examples just type inside the coinClpFBA folder:

```
make all
```

If the program compiled succesfully you can run the two binary examples by typing
```
./bin/run_fba data/iSIM.xml
```
or
```
./bin/run_fba data/iSIM.xml
```
To sbml_parse binary just read the SBML extract the info and print the data used to create the LP model:
```
./bin/parse_sbml data/Ecoli_core.xml
R_ACALD c=0 (-1000,1000)
  1 M_coa_c + 1 M_acald_c + 1 M_nad_c <==> 1 M_accoa_c + 1 M_h_c + 1 M_nadh_c
R_ACALDt c=0 (-1000,1000)
  1 M_acald_e <==> 1 M_acald_c
R_ACKr c=0 (-1000,1000)
  1 M_atp_c + 1 M_ac_c <==> 1 M_actp_c + 1 M_adp_c
...
```


Insted of the mode iSIM.xml you can also try Ecoli_core.xml.
However fo some strange reason I do not understand yet, FBA fails on bigger models such ass:
* iJO1366.xml
* Recon2.2.1_RPMI_trimed_gene_symbols.xml

## macOS

libSBML: (since the recommended way to build is to use [cmake](https://cmake.org/), that's what we do)
```
cd libSBML-5.16.0-Source
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_C_COMPILER=gcc-8 ..   # Use whatever OpenMP-enabled C++ compiler as used for PhysiCell
make
```
this should create `src/libsbml.dylib`
