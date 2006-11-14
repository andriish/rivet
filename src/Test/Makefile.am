noinst_PROGRAMS = testVector testLogging testHepPDT

MYLDFLAGS  = -L$(top_builddir)/src/

testVector_SOURCES   = testVector.cc
testVector_LDFLAGS   = $(MYLDFLAGS)
testVector_LDADD     = -lRivetTools

testLogging_SOURCES  = testLogging.cc
testLogging_LDFLAGS  = $(MYLDFLAGS) -L$(LOG4CPPLIBPATH)
testLogging_LDADD    = -l$(LOG4CPPLIBNAME)

testHepPDT_SOURCES   = testHepPDT.cc
testHepPDT_LDFLAGS   = $(MYLDFLAGS) -L$(LOG4CPPLIBPATH) -L$(HEPIDLIBPATH)
testHepPDT_LDADD     = -l$(LOG4CPPLIBNAME) -l$(HEPPIDLIBNAME)
testHepPDT_CPPFLAGS  = $(AM_CPPFLAGS) -I$(HEPPIDINCPATH)
