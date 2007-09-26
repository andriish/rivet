noinst_LTLIBRARIES = libRivetMath.la libRivetTools.la

libRivetMath_la_SOURCES  = 

libRivetTools_la_SOURCES  = Logging.cc RivetAIDA.cc Commandline.cc RivetPaths.cc \
                            TinyXML/tinyxml.cpp TinyXML/tinyxmlerror.cpp TinyXML/tinyxmlparser.cpp
libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) \
  -I$(top_srcdir)/include/TinyXML \
  -DLWH_USING_AIDA \
  -DRIVETDATADIR=\"$(datadir)\" \
  -DRIVETLIBDIR=\"$(libdir)\"
libRivetTools_la_LIBADD   = libRivetMath.la
if WITH_ROOT
  libRivetTools_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif
