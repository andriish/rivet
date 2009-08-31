lib_LTLIBRARIES = RivetRHICAnalyses.la

RivetRHICAnalyses_la_SOURCES = \
    STAR_2006_S6870392.cc \
    STAR_2008_S7993412.cc
    #STAR_2009_UE_HELEN.cc

## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
