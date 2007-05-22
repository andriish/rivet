noinst_LTLIBRARIES = libRivetMath.la libRivetTools.la #libRivetExtTools.la

libRivetMath_la_SOURCES  = 

#libRivetExtTools_la_SOURCES = 

libRivetTools_la_SOURCES  = Logging.cc RivetAIDA.cc Commandline.cc RivetPaths.cc \
                            TinyXML/tinyxml.cpp TinyXML/tinyxmlerror.cpp TinyXML/tinyxmlparser.cpp
libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) -DLWH_USING_AIDA -I$(top_srcdir)/include/TinyXML
libRivetTools_la_LIBADD   = libRivetMath.la


if WITH_ROOT
  libRivetTools_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS) -fpermissive
endif
