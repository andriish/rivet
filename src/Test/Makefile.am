noinst_PROGRAMS = testVector testLogging

MYLDFLAGS  = -L$(top_builddir)/src/

testVector_SOURCES  = testVector.cc
testVector_LDFLAGS  = $(MYLDFLAGS)
testVector_LDADD    = -lRivet

testLogging_SOURCES  = testLogging.cc
testLogging_LDFLAGS  = $(MYLDFLAGS) -L$(LOG4CPPLIBPATH)
testLogging_LDADD    = -l$(LOG4CPPLIBNAME)
