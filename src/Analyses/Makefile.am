lib_LTLIBRARIES = libRivetAnalysisStd.la

libRivetAnalysisStd_la_SOURCES = StdAnalyses.cc \
    TestAnalysis.cc \
    ExampleTree.cc  \
    HepEx9506012.cc \
    HepEx0112029.cc \
    HepEx0409040.cc \
    PL273B181.cc    \
    PRD65092002.cc  \
    ZPhys73C11.cc

libRivetAnalysisStd_la_CPPFLAGS = $(AM_CPPFLAGS)
if WITH_ROOT
  libRivetAnalysisStd_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif
