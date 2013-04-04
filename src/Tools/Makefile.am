noinst_LTLIBRARIES = libRivetTools.la

libRivetTools_la_SOURCES = \
  Logging.cc \
  RivetYODA.cc \
  RivetPaths.cc \
  binreloc.c \
  Utils.cc \
  ParticleIdUtils.cc \
  mt2_bisect.cc \
  RivetMT2.cc \
  BinnedHistogram.cc

dist_noinst_HEADERS = binreloc.h mt2_bisect.hh

libRivetTools_la_CPPFLAGS = \
  $(AM_CPPFLAGS) \
  -DENABLE_BINRELOC \
  -DDEFAULTDATADIR=\"$(datadir)\" \
  -DDEFAULTLIBDIR=\"$(libdir)\" \
  -I$(YAML_CPPINCDIR)

libRivetTools_la_LDFLAGS = $(AM_LDFLAGS) -L$(YAML_CPPLIBPATH)
libRivetTools_la_LIBADD = -lyaml-cpp
