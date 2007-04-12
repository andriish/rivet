noinst_LTLIBRARIES = libRivetAnalysis.la

libRivetAnalysis_la_SOURCES = \
	Analysis.cc  \
    TestAnalysis.cc \
    PL273B181.cc    \
	HepEx9506012.cc \
	HepEx0112029.cc \
	HepEx0409040.cc

## TODO We don't want this extra include polluting the include tree: need to remove it!
libRivetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) -I../../include/Rivet/Tools/D0RunIIcone
