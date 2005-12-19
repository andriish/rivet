LIBHERWIG = $(HOME)/local/lib/libherwig65.a
INCHERWIG = $(HOME)/local/include/herwig6507
LIBHEPMC = $(HOME)/local/lib/libHepMC.so
INCHEPMC = $(HOME)/local/include/
LIBCLHEP = /usr/local/lib/libCLHEP.so
INCCLHEP = /usr/local/include/
INCCERNLIB = $(HOME)
LIBCERNLIB = /usr/local/lib
#INCPDFLIB = $(HOME)/local/include/pdflib/
#LIBPDFLIB = $(HOME)/local/lib/libpdflib.a
INCPDFLIB = $(HOME)
LIBPDFLIB = $(HOME)/local/lib/libLHAPDF.a

DEBUG = -g
CC = g++
F77 = g77

CCHerwig: Makefile CCHerwig.o Herwig.o StringParse.o InitHerwigCommonBlocks_65.o herwig_common_block_address_65.o \
          HerwigInterface_65.o Address_i.o Address_f.o Address_c.o Address_cx.o upveto_hw.o topdec.o upinit_hw.o 
	$(CC) $(DEBUG) -I$(INCHEPMC) -I$(INCCLHEP) -I$(INCPDFLIB) \
	CCHerwig.cc Herwig.o StringParse.o InitHerwigCommonBlocks_65.o HerwigInterface_65.o \
	herwig_common_block_address_65.o \
	Address_i.o Address_f.o Address_c.o Address_cx.o upveto_hw.o upinit_hw.o topdec.o \
	$(LIBHERWIG) $(LIBCLHEP) $(LIBHEPMC) $(LIBPDFLIB) -lg2c \
	-o CCHerwig

CCHerwig.o: CCHerwig.cc
	$(CC) -c $(DEBUG) -I$(INCHEPMC) -I$(INCCLHEP) CCHerwig.cc

Herwig.o: Herwig.cxx
	$(CC) -c $(DEBUG) -I$(INCHEPMC) -I$(INCCLHEP) Herwig.cxx

InitHerwigCommonBlocks_65.o: InitHerwigCommonBlocks_65.cxx
	$(CC) -c $(DEBUG) InitHerwigCommonBlocks_65.cxx

StringParse.o: StringParse.cxx StringParse.h
	$(CC) -c $(DEBUG) StringParse.cxx

Address_i.o: Address_i.cxx
	$(CC) -c $(DEBUG) Address_i.cxx

Address_f.o: Address_f.cxx
	$(CC) -c $(DEBUG) Address_f.cxx

Address_c.o: Address_c.cxx
	$(CC) -c $(DEBUG) Address_c.cxx

Address_cx.o: Address_cx.cxx
	$(CC) -c $(DEBUG) Address_cx.cxx

upinit_hw.o: upinit_hw.F
	$(F77) -c $(DEBUG) -I$(INCHERWIG) upinit_hw.F

upveto_hw.o: upveto_hw.F
	$(F77) -c $(DEBUG) -I$(INCHERWIG) upveto_hw.F

topdec.o: topdec.F
	$(F77) -c $(DEBUG) -I$(INCHERWIG) topdec.F

HerwigInterface_65.o: HerwigInterface_65.F
	$(F77) -c $(DEBUG) -I$(INCHERWIG) HerwigInterface_65.F

herwig_common_block_address_65.o: herwig_common_block_address_65.F
	$(F77) -c $(DEBUG) -I$(INCHERWIG) herwig_common_block_address_65.F

clean:
	@rm -f *.o a.out
