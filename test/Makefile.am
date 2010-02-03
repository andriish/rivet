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

export RIVET_REF_PATH=$(top_builddir)/data/refdata
export RIVET_DATA_PATH=$(top_builddir)/data/refdata
export LD_LIBRARY_PATH=$(top_builddir)/src:$(top_builddir)/src/.libs:$(top_builddir)/src/Analyses/.libs:$(HEPMCLIBPATH):$(FASTJETLIBPATH)
export PYTHONPATH=$(shell ls -d $(top_builddir)/pyext/build/lib.*)
export PATH=$(top_builddir)/bin:$(shell echo $$PATH)

TESTS = testMatVec testBoost testCmp testApi testCmdLine.sh

EXTRA_DIST = testApi.hepmc
