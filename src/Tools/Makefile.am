noinst_LTLIBRARIES = libRivetTools.la

libRivetTools_la_SOURCES  = Logging.cc RivetAIDA.cc RivetPaths.cc binreloc.c Utils.cc ParticleIDMethods.cc \
                            TinyXML/tinyxml.cpp TinyXML/tinyxmlerror.cpp TinyXML/tinyxmlparser.cpp
libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) \
  -I$(top_srcdir)/include/TinyXML \
  -DLWH_USING_AIDA -DENABLE_BINRELOC \
  -DDEFAULTDATADIR=\"$(datadir)\" \
  -DDEFAULTLIBDIR=\"$(libdir)\"

if WITH_ROOT
  libRivetTools_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif
