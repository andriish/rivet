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

export RIVET_REF_PATH=$(datadir)/Rivet
export RIVET_DATA_PATH=$(datadir)/Rivet

TESTS = testMatVec testBoost testCmp testApi

EXTRA_DIST = testApi.hepmc
