noinst_PROGRAMS = testLogging testRivetHandler testHepPDT testLWH testTinyXML 

MYLDFLAGS  = -L$(top_builddir)/src -L$(top_builddir)/src/Tools

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
testTinyXML_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)/include

check-local:
	./testLogging
	./testRivetHandler
	./testHepPDT
	./testLWH && test -e test.aida.xml && rm -f test.aida.xml
#	./testTinyXML ../../data/4781.hepml.xml
