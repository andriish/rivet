lib_LTLIBRARIES = RivetD0Analyses.la

RivetD0Analyses_la_SOURCES = \
    D0_1996_S3214044.cc \
    D0_1996_S3324664.cc \
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

## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
