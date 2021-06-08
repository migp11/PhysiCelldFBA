After unpacking coin-or run the following comand:

# For some extrange reason in "coin-or/include/coin-or/ClpSimplex.hpp" and "coin-or/include/coin-or/ClpModel.hpp"
# the include to glpk.h has a wrong path

ln ext/coin-or/include/coin-or/glpk/glpk.h ext/coin-or/include/coin-or/

after compiling we need to export LD

---------------------------------------------------

# Install Anaconda
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh

chmod +x Anaconda3-2021.05-Linux-x86_64.sh 
./Anaconda3-2021.05-Linux-x86_64.sh 

export PATH=$PATH:./anaconda3/bin/
conda init bash

# create a conda environment
conda create --name pc-devel-env

# activate the conda environment
conda activate pc-devel-env

# installing libsbml
conda install -c vincent-noel libsbml-plus-packages

# installing coin-cbc MILP solver
conda install -c conda-forge coincbc


cd <path/to/addon>

make clean

g++ -Wall -std=c++11 -pipe -L/home/mponce/anaconda3/envs/pc-devel-env/lib -I/home/mponce/anaconda3/envs/pc-devel-env/include/ -o ./bin/testSBML ./test/testSBML.cpp -lsbml

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/anaconda3/envs/pc-devel-env/lib/

./bin/testSBML test/data/Ecoli_core.xml 
------------------------------------------------------------------------- 


Links:


libSBML-5.19.0-Linux-x64-binaries-ubuntu.tar.gz
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/stable/Linux/64-bit/libSBML-5.19.0-Linux-x64-binaries-ubuntu.tar.gz/download

libSBML-5.19.0-Linux-x64-binaries-centos.tar.gz
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/stable/Linux/64-bit/libSBML-5.19.0-Linux-x64-binaries-centos.tar.gz/download

libSBML-5.19.0-Linux-x86-binaries-ubuntu.tar.gz
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/stable/Linux/32-bit/libSBML-5.19.0-Linux-x86-binaries-ubuntu.tar.gz/download

libSBML-5.19.0-Linux-x86-binaries-centos.tar.gz
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/stable/Linux/32-bit/libSBML-5.19.0-Linux-x86-binaries-centos.tar.gz/download


libSBML-5.19.0-packages-Linux-x86-ubuntu-binaries.tar.gz
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/experimental/binaries/Linux/libSBML-5.19.0-packages-Linux-x86-ubuntu-binaries.tar.gz/download


libSBML-5.19.0-Darwin.tar.gz
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/experimental/binaries/Mac%20OS%20X/libSBML-5.19.0-Darwin.tar.gz/download

libSBML-5.19.0-win64.zip
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/experimental/binaries/Windows/libSBML-5.19.0-win64.zip/download

libSBML-5.19.0-win32.zip
https://sourceforge.net/projects/sbml/files/libsbml/5.19.0/experimental/binaries/Windows/libSBML-5.19.0-win32.zip/download
