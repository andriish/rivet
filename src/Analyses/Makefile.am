## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
lib_LTLIBRARIES =


## Unvalidated analyses: add new standard analyses here, and only
## move them into the collider-specific standard plugin libraries
## once they have been finished and checked. The --enable-unvalidated
## flag is needed to install the RivetUnvalidated library.
if ENABLE_UNVALIDATED
lib_LTLIBRARIES += RivetUnvalidated.la
RivetUnvalidated_la_SOURCES = \
    CDF_1991_S2313472.cc \
    CDF_1996_S3108457.cc \
    CDF_1996_S3349578.cc \
    CDF_1996_S3418421.cc \
    CDF_2000_S4266730.cc \
    CDF_2001_S4517016.cc \
    CDF_2001_S4563131.cc \
    CDF_1994_S2952106.cc \
    CDF_2006_S6653332.cc \
    CDF_2007_S7057202.cc \
    CDF_2008_S7541902.cc \
    CDF_2008_S7782535.cc \
    CDF_2008_S8095620.cc \
    D0_1996_S3214044.cc \
    D0_1996_S3324664.cc \
    D0_1998_S3711838.cc \
    DELPHI_2003_WUD_03_11.cc \
    E735_1998_S3905616.cc \
    H1_1995_S3167097.cc \
    MC_DIJET.cc \
    MC_PHOTONJETUE.cc \
    MC_TTBAR.cc \
    SFM_1984_S1178091.cc \
    STAR_2008_S7993412.cc \
    STAR_2009_UE_HELEN.cc \
    UA1_1990_S2044935.cc \
    UA5_1988_S1867512.cc \
    ZEUS_2001_S4815815.cc
endif


lib_LTLIBRARIES += RivetCDFAnalyses.la
RivetCDFAnalyses_la_SOURCES = \
    CDF_1988_S1865951.cc \
    CDF_1990_S2089246.cc \
    CDF_1998_S3618439.cc \
    CDF_2000_S4155203.cc \
    CDF_2001_S4751469.cc \
    CDF_2002_S4796047.cc \
    CDF_2004_S5839831.cc \
    CDF_2005_S6080774.cc \
    CDF_2005_S6217184.cc \
    CDF_2006_S6450792.cc \
    CDF_2008_LEADINGJETS.cc \
    CDF_2008_NOTE_9351.cc \
    CDF_2008_S7540469.cc \
    CDF_2008_S7828950.cc \
    CDF_2008_S8093652.cc \
    CDF_2009_S8233977.cc \
    CDF_2009_S8383952.cc \
    CDF_2009_S8436959.cc


lib_LTLIBRARIES += RivetD0Analyses.la
RivetD0Analyses_la_SOURCES = \
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
    D0_2009_S8202443.cc \
    D0_2009_S8320160.cc \
    D0_2009_S8349509.cc


lib_LTLIBRARIES += RivetHERAAnalyses.la
RivetHERAAnalyses_la_SOURCES = \
    H1_1994_S2919893.cc \
    H1_2000_S4129130.cc


lib_LTLIBRARIES += RivetLEPAnalyses.la
RivetLEPAnalyses_la_SOURCES = \
    ALEPH_1991_S2435284.cc \
    ALEPH_1996_S3486095.cc \
    ALEPH_1996_S3196992.cc \
    DELPHI_1995_S3137023.cc \
    DELPHI_1996_S3430090.cc \
    DELPHI_2002_069_CONF_603.cc \
    OPAL_1998_S3780481.cc \
    OPAL_2004_S6132243.cc


lib_LTLIBRARIES += RivetRHICAnalyses.la
RivetRHICAnalyses_la_SOURCES = \
    STAR_2006_S6500200.cc \
    STAR_2006_S6860818.cc \
    STAR_2006_S6870392.cc


lib_LTLIBRARIES += RivetSPSAnalyses.la
RivetSPSAnalyses_la_SOURCES = \
    UA5_1982_S875503.cc \
    UA5_1986_S1583476.cc \
    UA5_1989_S1926373.cc


lib_LTLIBRARIES += RivetMiscAnalyses.la
RivetMiscAnalyses_la_SOURCES = \
    PDG_Hadron_Multiplicities.cc \
    PDG_Hadron_Multiplicities_Ratios.cc \
    JADE_OPAL_2000_S4300807.cc \
    BELLE_2006_S6265367.cc


lib_LTLIBRARIES += RivetMCAnalyses.la
RivetMCAnalyses_la_SOURCES = \
    ExampleAnalysis.cc \
    MC_JetAnalysis.cc \
    MC_PHOTONJETS.cc \
    MC_ZJETS.cc \
    MC_LEADINGJETS.cc \
    MC_DIPHOTON.cc \
    MC_WANALYSIS.cc \
    MC_ZANALYSIS.cc \
    MC_SUSY.cc
