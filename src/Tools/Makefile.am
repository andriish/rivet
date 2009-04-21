noinst_LTLIBRARIES = libRivetTools.la

libRivetTools_la_SOURCES  = Logging.cc RivetAIDA.cc RivetPaths.cc binreloc.c Utils.cc ParticleIDMethods.cc \
                            TinyXML/tinyxml.cpp TinyXML/tinyxmlerror.cpp TinyXML/tinyxmlparser.cpp \
                            BinnedHistogram.cc

dist_noinst_HEADERS = binreloc.h

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) \
  -I$(top_srcdir)/include/TinyXML \
  -DENABLE_BINRELOC \
  -DDEFAULTDATADIR=\"$(datadir)\" \
  -DDEFAULTLIBDIR=\"$(libdir)\"
