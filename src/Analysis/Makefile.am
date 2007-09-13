noinst_LTLIBRARIES = libRivetAnalysis.la
lib_LTLIBRARIES = libRivetAnalysisStd.la

libRivetAnalysis_la_SOURCES = Analysis.cc AnalysisLoader.cc
libRivetAnalysisStd_la_SOURCES = StdAnalyses.cc \
    TestAnalysis.cc \
    ExampleTree.cc  \
    HepEx9506012.cc \
    HepEx0112029.cc \
    HepEx0409040.cc \
    PL273B181.cc    \
    PRD65092002.cc

if WITH_ROOT
  libRivetAnalysisStd_la_CPPFLAGS = $(AM_CPPFLAGS) -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
else 
  libRivetAnalysisStd_la_CPPFLAGS = $(AM_CPPFLAGS)
endif
