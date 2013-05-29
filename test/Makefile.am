check_PROGRAMS = testMath testMatVec testBoost testCmp testApi testNaN

AM_LDFLAGS = -L$(top_srcdir)/src -L$(prefix)/lib $(BOOST_LDFLAGS) $(YAMLCPP_LDFLAGS) -L$(YODALIBPATH)
LIBS = -lm $(GSL_LDFLAGS) -lYODA

testBoost_SOURCES = testBoost.cc
testMath_SOURCES = testMath.cc
testMatVec_SOURCES = testMatVec.cc

testCmp_SOURCES = testCmp.cc
testCmp_LDADD = $(top_builddir)/src/libRivet.la $(HEPMCLDFLAGS) $(HEPMCLDLIBS)
testApi_SOURCES = testApi.cc
testApi_LDADD = $(top_builddir)/src/libRivet.la $(HEPMCLDFLAGS) $(HEPMCLDLIBS)
testNaN_SOURCES = testNaN.cc
testNaN_LDADD = $(top_builddir)/src/libRivet.la $(HEPMCLDFLAGS) $(HEPMCLDLIBS)

TESTS_ENVIRONMENT = \
  RIVET_REF_PATH=$(top_builddir)/data/refdata \
  RIVET_INFO_PATH=$(top_builddir)/data/anainfo \
  RIVET_ANALYSIS_PATH=$(top_builddir)/src/Analyses/.libs \
  LD_LIBRARY_PATH=$(top_builddir)/src:$(top_builddir)/src/.libs:$(HEPMCLIBPATH):$(FASTJETLIBPATH):$(YAML_CPPLIBPATH):$(LD_LIBRARY_PATH) \
  DYLD_LIBRARY_PATH=$(top_builddir)/src:$(top_builddir)/src/.libs:$(HEPMCLIBPATH):$(FASTJETLIBPATH):$(YAML_CPPLIBPATH):$(DYLD_LIBRARY_PATH) \
  PYTHON_BUILD_DIR=$(top_builddir)/pyext/build \
  PATH=$(top_builddir)/bin:$(PATH)

TESTS = testMath testMatVec testBoost testCmp testApi testNaN testCmdLine.sh

EXTRA_DIST = testApi.hepmc testCmdLine.sh

clean-local:
	@rm -f out.aida log a.out fifo.hepmc file2.hepmc
