noinst_PROGRAMS = testVector testLogging testRivetHandler testHepPDT testLWH

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
#-l$(HEPPIDLIBNAME) -lRivetTools
testHepPDT_CPPFLAGS  = $(AM_CPPFLAGS) -I$(HEPPIDINCPATH)

testLWH_SOURCES  = testLWH.cc
#testLogging_LDFLAGS  = $(MYLDFLAGS)
#testLogging_LDADD    = -lRivet


check-local:
	./testVector
	./testLogging
	./testRivetHandler
	./testHepPDT
