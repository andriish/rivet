#! /usr/bin/env bash

LIBHERWIG=$HOME/local/lib/libherwig65.a
INCHERWIG=$HOME/local/include/herwig6507
LIBHEPMC=$HOME/local/lib/libHepMC.so
INCHEPMC=$HOME/local/include/
LIBCLHEP=/usr/local/lib/libCLHEP.so
INCCLHEP=/usr/local/include/

g++ -I$INCHEPMC -I$INCCLHEP Herwig.cxx 2>&1 
