lib_LTLIBRARIES = libRivetAnalysisStd.la

libRivetAnalysisStd_la_SOURCES = StdAnalyses.cc \
    TestAnalysis.cc \
    ExampleTree.cc  \
    S3167097.cc \
    HepEx0107012.cc \
    S4815815.cc \
    S5992206.cc \
    S6217184.cc \
    S7057202.cc \
    S2435284.cc    \
    S4751469.cc  \
    S3430090.cc   \
    HepEx0605099.cc 
#    MyAnalysis.cc	

libRivetAnalysisStd_la_CPPFLAGS = $(AM_CPPFLAGS)
if WITH_ROOT
  libRivetAnalysisStd_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif
