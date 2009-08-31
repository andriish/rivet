lib_LTLIBRARIES = RivetMCAnalyses.la

RivetMCAnalyses_la_SOURCES = \
    MC_JetAnalysis.cc \
    MC_TVT1960_PHOTONJETS.cc \
    MC_TVT1960_ZJETS.cc \
    MC_LHC_LEADINGJETS.cc \
    MC_LHC_DIJET.cc \
    MC_LHC_ZANALYSIS.cc

## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
