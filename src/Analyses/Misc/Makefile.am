lib_LTLIBRARIES = RivetMiscAnalyses.la

RivetMiscAnalyses_la_SOURCES = \
    PDG_Hadron_Multiplicities.cc \
    PDG_Hadron_Multiplicities_Ratios.cc \
    SFM_1984_S1178091.cc \
    E735_1998_S3905616.cc \
    JADE_OPAL_2000_S4300807.cc

## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
