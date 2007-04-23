noinst_LTLIBRARIES = libRivetMath.la libRivetExtTools.la libRivetTools.la

libRivetMath_la_SOURCES  = 

libRivetExtTools_la_SOURCES = 

libRivetTools_la_SOURCES  = Logging.cc RivetAIDA.cc Commandline.cc RivetPaths.cc
libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) -DLWH_USING_AIDA -I$(top_srcdir)/include/TinyXML
libRivetTools_la_LIBADD   = \
    libRivetMath.la libRivetXML.la libRivetExtTools.la \
	TinyXML/tinyxml.cpp TinyXML/tinyxmlerror.cpp TinyXML/tinyxmlparser.cpp
