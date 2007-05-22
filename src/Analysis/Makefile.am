noinst_LTLIBRARIES = libRivetAnalysis.la

libRivetAnalysis_la_SOURCES = \
	Analysis.cc     \
	TestAnalysis.cc \
        ExampleTree.cc \
	HepEx9506012.cc \
	HepEx0112029.cc \
	HepEx0409040.cc \
        PL273B181.cc    \
        PRD65092002.cc  

if WITH_ROOT
  libRivetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
else 
  libRivetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
endif

