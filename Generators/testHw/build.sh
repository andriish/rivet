LIBHERWIG=$(HOME)/local/lib/libherwig65.a
INCHERWIG=$(HOME)/local/include/herwig6507
LIBHEPMC=$(HOME)/local/lib/libHepMC.so
INCHEPMC=$(HOME)/local/include/
LIBCLHEP=/usr/local/lib/libCLHEP.so
INCCLHEP=/usr/local/include/

all:
	g++ -I$(INCHEPMC) -I$(INCCLHEP) $(LIBHERWIG) $(LIBCLHEP) Herwig.cxx 2>&1
