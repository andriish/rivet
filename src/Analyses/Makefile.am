lib_LTLIBRARIES = libRivetAnalysisStd.la

libRivetAnalysisStd_la_SOURCES = StdAnalyses.cc \
    ExampleAnalysis.cc \
    ExampleTree.cc  \
    S2435284.cc \
    S3167097.cc \
    S3430090.cc \
    S4674421.cc \
    S4751469.cc \
    S4815815.cc \
    S5992206.cc \
    S6217184.cc \
    S6653332.cc \
    S7057202.cc

libRivetAnalysisStd_la_CPPFLAGS = $(AM_CPPFLAGS)
if WITH_ROOT
  libRivetAnalysisStd_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif
