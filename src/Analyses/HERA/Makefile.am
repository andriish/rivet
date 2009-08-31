lib_LTLIBRARIES = RivetHERAAnalyses.la

RivetHERAAnalyses_la_SOURCES = \
    H1_1994_S2919893.cc \
    H1_1995_S3167097.cc \
    H1_2000_S4129130.cc \
    ZEUS_2001_S4815815.cc

## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
