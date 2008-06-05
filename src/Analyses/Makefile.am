lib_LTLIBRARIES = libRivetAnalysisStd.la

libRivetAnalysisStd_la_SOURCES = StdAnalyses.cc \
    ExampleAnalysis.cc \
    ExampleTree.cc  \
    ALEPH_1991_S2435284.cc \
    CDF_1994_S2952106.cc \
    CDF_2001_S4751469.cc \
    CDF_2004_S5839831.cc \
    CDF_2005_S6217184.cc \
    CDF_2006_S6653332.cc \
    CDF_2007_S7057202.cc \
    CDF_2008_S7541902.cc \
    D0_2001_S4674421.cc \
    D0_2004_S5992206.cc \
    DELPHI_1995_S3137023.cc \
    DELPHI_1996_S3430090.cc \
    DELPHI_2002_069_CONF_603.cc \
    H1_1995_S3167097.cc \
    OPAL_1998_S3780481.cc \
    ZEUS_2001_S4815815.cc \
    PDG_Hadron_Multiplicities.cc
    #OPAL_2004_S6132243.cc 

libRivetAnalysisStd_la_LDFLAGS = $(AM_LDFLAGS) -L$(FASTJETLIBPATH)
libRivetAnalysisStd_la_LIBADD = -lCDFConesPlugin -lSISConePlugin -lsiscone
libRivetAnalysisStd_la_CPPFLAGS = $(AM_CPPFLAGS)

if WITH_ROOT
  libRivetAnalysisStd_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif

if ENABLE_JADE
  libRivetAnalysisStd_la_CPPFLAGS += -DHAVE_JADE
  libRivetAnalysisStd_la_LIBADD += -lJadePlugin
endif
