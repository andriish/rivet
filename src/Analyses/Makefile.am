lib_LTLIBRARIES = libRivetAnalysisStd.la

libRivetAnalysisStd_la_SOURCES = StdAnalyses.cc \
    ExampleAnalysis.cc \
    ExampleTree.cc  \
    ALEPH_1991_S2435284.cc \
    ALEPH_1996_S3486095.cc \
    CDF_1994_S2952106.cc \
    CDF_2000_S4155203.cc \
    CDF_2001_S4751469.cc \
    CDF_2002_S4796047.cc \
    CDF_2004_S5839831.cc \
    CDF_2005_S6217184.cc \
    CDF_2006_S6653332.cc \
    CDF_2007_S7057202.cc \
    CDF_2008_S7541902.cc \
    CDF_2008_NOTE_9337.cc \
    CDF_2008_NOTE_9351.cc \
    CDF_2008_LEADINGJETS.cc \
    CDF_2008_S7540469.cc \
    CDF_2008_S7782535.cc \
    CDF_2008_S7828950.cc \
    D0_2001_S4674421.cc \
    D0_2004_S5992206.cc \
    D0_2006_S6438750.cc \
	D0_2007_S7075677.cc \
    D0_2008_S6879055.cc \
    D0_2008_S7554427.cc \
    D0_2008_S7662670.cc \
    D0_2008_S7719523.cc \
    D0_2008_S7837160.cc \
    D0_2008_S7863608.cc \
    DELPHI_1995_S3137023.cc \
    DELPHI_1996_S3430090.cc \
    DELPHI_2002_069_CONF_603.cc \
    DELPHI_2003_WUD_03_11.cc \
    H1_1994_S2919893.cc \
    H1_1995_S3167097.cc \
    H1_2000_S4129130.cc \
    OPAL_1998_S3780481.cc \
    ZEUS_2001_S4815815.cc \
    PDG_Hadron_Multiplicities.cc \
    PDG_Hadron_Multiplicities_Ratios.cc 
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
