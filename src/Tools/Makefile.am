noinst_LTLIBRARIES = libRivetMath.la libRivetXML.la libRivetTools.la

#libRivetMath_la_SOURCES  = 

libRivetXML_la_SOURCES   = TinyXML/tinyxml.cpp TinyXML/tinyxmlerror.cpp TinyXML/tinyxmlparser.cpp
libRivetXML_la_CPPFLAGS  = -I$(top_srcdir)/include/TinyXML

libRivetTools_la_SOURCES  = Logging.cc RivetAIDA.cc
libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) -DLWH_USING_AIDA
libRivetTools_la_LDFLAGS  = -L$(LOG4CPPLIBPATH) -R$(LOG4CPPLIBPATH)
libRivetTools_la_LIBADD   = libRivetMath.la libRivetXML.la -l$(LOG4CPPLIBNAME)
