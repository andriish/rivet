noinst_PROGRAMS = testVector testLogging testHepPDT

MYLDFLAGS  = -L$(top_builddir)/src -L$(top_builddir)/src/Tools

testVector_SOURCES   = testVector.cc
testVector_LDFLAGS   = $(MYLDFLAGS)
testVector_LDADD     = -lRivet

testLogging_SOURCES  = testLogging.cc
testLogging_LDFLAGS  = $(MYLDFLAGS)
testLogging_LDADD    = -lRivet

testHepPDT_SOURCES   = testHepPDT.cc
testHepPDT_LDFLAGS   = $(MYLDFLAGS) -L$(HEPIDLIBPATH)
testHepPDT_LDADD     = -lRivet 
#-l$(HEPPIDLIBNAME) -lRivetTools
testHepPDT_CPPFLAGS  = $(AM_CPPFLAGS) -I$(HEPPIDINCPATH)
