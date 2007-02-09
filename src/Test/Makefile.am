noinst_PROGRAMS = testVector testLogging testRivetHandler testHepPDT testLWH testTinyXML 

MYLDFLAGS  = -L$(top_builddir)/src -L$(top_builddir)/src/Tools

testVector_SOURCES   = testVector.cc
testVector_LDFLAGS   = $(MYLDFLAGS)
testVector_LDADD     = -lRivet -lHepPDT

testLogging_SOURCES  = testLogging.cc
testLogging_LDFLAGS  = $(MYLDFLAGS)
testLogging_LDADD    = -lRivet

testRivetHandler_SOURCES  = testRivetHandler.cc
testRivetHandler_LDFLAGS  = $(MYLDFLAGS)
testRivetHandler_LDADD    = -lRivet

testHepPDT_SOURCES   = testHepPDT.cc
testHepPDT_LDFLAGS   = $(MYLDFLAGS) -L$(HEPIDLIBPATH)
testHepPDT_LDADD     = -lRivet 
testHepPDT_CPPFLAGS  = $(AM_CPPFLAGS) -I$(HEPPIDINCPATH)

testLWH_SOURCES  = testLWH.cc
testLWH_CPPFLAGS = $(AM_CPPFLAGS) -DLWH_USING_AIDA

testTinyXML_SOURCES  = testTinyXML.cc
testTinyXML_LDFLAGS  = $(MYLDFLAGS)
testTinyXML_LDADD    = -lRivetXML
testTinyXML_CPPFLAGS = $(AM_CPPFLAGS)

if WITH_XERCES_C
noinst_PROGRAMS += testXerces 
testXerces_SOURCES  = testXerces.cc
testXerces_CPPFLAGS = $(AM_CPPFLAGS) -I$(XERCES_CINCPATH)
testXerces_LDFLAGS = $(AM_LDFLAGS) -L$(XERCES_CLIBPATH)
testXerces_LDADD = -l$(XERCES_CLIBNAME)
endif

check-local:
	./testVector
	./testLogging
	./testRivetHandler
	./testHepPDT
	./testLWH
	./testTinyXML
