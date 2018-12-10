noinst_LTLIBRARIES = libRivetTools.la

libRivetTools_la_SOURCES = \
  BinnedHistogram.cc \
  Cuts.cc \
  JetUtils.cc \
  Random.cc \
  Logging.cc \
  ParticleUtils.cc \
  ParticleName.cc \
  RivetYODA.cc \
  RivetMT2.cc \
  RivetPaths.cc \
  Utils.cc \
  binreloc.c
#  mt2_bisect.cc

if ENABLE_HEPMC_3
libRivetTools_la_SOURCES += RivetHepMC_3.cc
else
libRivetTools_la_SOURCES += RivetHepMC_2.cc
endif

dist_noinst_HEADERS = binreloc.h lester_mt2_bisect.hh mt2_bisect.hh

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) -DENABLE_BINRELOC -DDEFAULTDATADIR=\"$(datadir)\" -DDEFAULTLIBDIR=\"$(libdir)\"
