lib_LTLIBRARIES = RivetSPSAnalyses.la

RivetSPSAnalyses_la_SOURCES = \
    UA1_1990_S2044935.cc \
    UA5_1982_S875503.cc \
    UA5_1986_S1583476.cc \
    UA5_1988_S1867512.cc \
    UA5_1989_S1926373.cc

## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
