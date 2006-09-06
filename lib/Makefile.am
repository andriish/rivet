## Workaround for OS X, automake doesn't realize it's meant to use c++
LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
	$(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

pkglib_LTLIBRARIES = libRivet.la

# Version info should be updated if any interface or persistent I/O
# function is changed
libRivet_la_LDFLAGS = -version-info 1:0:0 -export-dynamic -L$(CLHEPLIBPATH)
libRivet_la_SOURCES =
libRivet_la_LIBADD = \
	$(top_builddir)/src/Analysis/libRivetAnalysis.la \
	$(top_builddir)/src/Analysis/libRivetAnalysisExamples.la \
	$(top_builddir)/src/Projections/libRivetProjections.la \
	-l$(CLHEPLIBNAME) -ldl 
