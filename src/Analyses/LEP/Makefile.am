lib_LTLIBRARIES = RivetLEPAnalyses.la

RivetLEPAnalyses_la_SOURCES = \
    ALEPH_1991_S2435284.cc \
    ALEPH_1996_S3486095.cc \
    DELPHI_1995_S3137023.cc \
    DELPHI_1996_S3430090.cc \
    DELPHI_2002_069_CONF_603.cc \
    DELPHI_2003_WUD_03_11.cc \
    OPAL_1998_S3780481.cc
    #OPAL_2004_S6132243.cc 

## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
