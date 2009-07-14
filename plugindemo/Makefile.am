## Build rule for the demo plugin analysis, MyAnalysis

if WITH_OSX
  pluginlib = libRivetMyAnalysis.dylib
  shared_flags =-undefined dynamic_lookup
else
  pluginlib = libRivetMyAnalysis.so
  shared_flags = -shared
endif

all:
	@echo "Default make rule does nothing: use 'make plugin'"

EXTRA_DIST = MyAnalysis.cc
$(pluginlib): MyAnalysis.cc
	$(CXX) -o $(pluginlib) $(shared_flags) -I$(top_srcdir)/include -I$(HEPMCINCPATH) MyAnalysis.cc

## Hook in the build rule for the plugin library demo
plugin: $(pluginlib)
	@true

## Test the plugin
EXTRA_DIST += testplugin
checkplugin: plugin
	RIVET_ANALYSIS_PATH=$(srcdir) $(srcdir)/testplugin
