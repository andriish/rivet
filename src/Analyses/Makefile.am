lib_LTLIBRARIES = libRivetAnalysisStd.la

libRivetAnalysisStd_la_SOURCES = StdAnalyses.cc \
    ExampleAnalysis.cc \
    ExampleTree.cc  \
    ALEPH_1991_S2435284.cc \
    CDF_2001_S4751469.cc \
    CDF_2005_S6217184.cc \
    CDF_2006_S6653332.cc \
    CDF_2007_S7057202.cc \
    D0_2001_S4674421.cc \
    D0_2004_S5992206.cc \
    DELPHI_1996_S3430090.cc \
    H1_1995_S3167097.cc \
    ZEUS_2001_S4815815.cc
    #OPAL_2004_S6132243.cc 

libRivetAnalysisStd_la_CPPFLAGS = $(AM_CPPFLAGS)
if WITH_ROOT
  libRivetAnalysisStd_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif
