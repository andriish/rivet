check_PROGRAMS = testMatVec testBoost testCmp testApi

testBoost_SOURCES = testBoost.cc
testBoost_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)

testMatVec_SOURCES = testMatVec.cc
testMatVec_LDADD   = -lm -lgsl -lgslcblas
testMatVec_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)

testCmp_SOURCES = testCmp.cc
testCmp_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)

testApi_SOURCES = testApi.cc
testApi_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)
testApi_LDADD = ../src/libRivet.la $(HEPMCLDFLAGS) $(HEPMCLDLIBS)

## TODO: Move these non-POSIX shell escapes into variable exports from configure.ac
TESTS_ENVIRONMENT = \
  RIVET_REF_PATH=$(top_builddir)/data/refdata \
  RIVET_INFO_PATH=$(top_builddir)/data/anainfo \
  RIVET_ANALYSIS_PATH=$(top_builddir)/src/Analyses/.libs \
  LD_LIBRARY_PATH=$(top_builddir)/src:$(top_builddir)/src/.libs:$(HEPMCLIBPATH):$(FASTJETLIBPATH) \
  PYTHONPATH=$(shell ls -d $(top_builddir)/pyext/build/lib.*) \
  PATH=$(top_builddir)/bin:$(shell echo $$PATH)

TESTS = testMatVec testBoost testCmp testApi testCmdLine.sh

EXTRA_DIST = testApi.hepmc testCmdLine.sh
