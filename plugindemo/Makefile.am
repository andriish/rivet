## Build rule for the demo plugin analysis, MyAnalysis
pluginlib = libRivetMyAnalysis$(RAW_LIB_EXTN)
$(pluginlib): MyAnalysis.cc
	$(CXX) -o $(pluginlib) -shared -I$(top_srcdir)/include -I$(HEPMCINCPATH) MyAnalysis.cc


## Hook in the build rule for the plugin library demo
all-local: $(pluginlib)
	@true


## Use rivet to test for the new plugin being loaded
if ENABLE_PYRIVET
EXTRA_DIST = testplugin
TESTS_ENVIRONMENT = RIVET_ANALYSIS_PATH=$(top_srcdir)/plugindemo
TESTS = testplugin
endif
