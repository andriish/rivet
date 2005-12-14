LIBHERWIG=$(HOME)/local/lib/libherwig65.a
INCHERWIG=$(HOME)/local/include/herwig6507
LIBHEPMC=$(HOME)/local/lib/libHepMC.so
INCHEPMC=$(HOME)/local/include/
LIBCLHEP=/usr/local/lib/libCLHEP.so
INCCLHEP=/usr/local/include/
INCCERNLIB=
LIBCERNLIB=/usr/local/lib

CCHerwig: StringParse.o InitHerwigCommonBlocks_65.o herwig_common_block_address_65.o HerwigInterface_65.o Address_i.o Address_f.o Address_c.o Address_cx.o Herwig.o
	g++ -I$(INCHEPMC) -I$(INCCLHEP) $(LIBHERWIG) $(LIBCLHEP) $(LIBHEPMC) -lg2c \
	Herwig.o StringParse.o InitHerwigCommonBlocks_65.o herwig_common_block_address_65.o HerwigInterface_65.o \
	Address_i.o Address_f.o Address_c.o Address_cx.o \
	CCHerwig.cc -o CCHerwig

Herwig.o:
	g++ -c -I$(INCHEPMC) -I$(INCCLHEP) Herwig.cxx

InitHerwigCommonBlocks_65.o:
	g++ -c InitHerwigCommonBlocks_65.cxx

StringParse.o: StringParse.cxx StringParse.h
	g++ -c StringParse.cxx

Address_i.o:
	g++ -c Address_i.cxx

Address_f.o:
	g++ -c Address_f.cxx

Address_c.o:
	g++ -c Address_c.cxx

Address_cx.o:
	g++ -c Address_cx.cxx

HerwigInterface_65.o:
	g77 -c -I$(INCHERWIG) HerwigInterface_65.F

herwig_common_block_address_65.o:
	g77 -c -I$(INCHERWIG) herwig_common_block_address_65.F

clean:
	@rm -f *.o a.out
