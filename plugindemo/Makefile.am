#noinst_LIBRARIES = libRivetMyAnalysis.la
#libRivetMyAnalysis_la_SOURCES = MyAnalysis.cc
#libRivetMyAnalysis_la_CPPFLAGS = -I$(top_srcdir)/include -I$(HEPMCINCPATH)

libRivetMyAnalysis.so: MyAnalysis.cc
	$(CXX) -o libRivetMyAnalysis.so -shared -I$(top_srcdir)/include -I$(HEPMCINCPATH) MyAnalysis.cc

if ENABLE_PYRIVET
EXTRA_DIST = testplugin
TESTS_ENVIRONMENT = RIVET_ANALYSIS_PATH=$(top_srcdir)/plugindemo
TESTS = testplugin
endif
