## Build rule for the demo plugin analysis, MyAnalysis

pluginlib = RivetMyAnalysis.so
if WITH_OSX
  shared_flags = -undefined dynamic_lookup -bundle 
else
  shared_flags = -shared -fPIC
endif

all:
	@echo "Default make rule does nothing: use 'make plugin'"

EXTRA_DIST = MyAnalysis.cc
$(pluginlib): MyAnalysis.cc
	$(CXX) -o $(pluginlib) $(shared_flags) -I$(top_srcdir)/include -I$(HEPMCINCPATH) MyAnalysis.cc

## Bundle the ROOT-based analysis in the Rivet tarball
EXTRA_DIST += RootAnalysis.cc

## Hook in the build rule for the plugin library demo
plugin: $(pluginlib)
	@true

## Test the plugin
EXTRA_DIST += testplugin
checkplugin: plugin
	RIVET_ANALYSIS_PATH=$(srcdir) $(srcdir)/testplugin
