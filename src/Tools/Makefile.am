noinst_LTLIBRARIES = libRivetMath.la libRivetEvent.la libRivetTools.la

libRivetMath_la_SOURCES  = Math/LorentzVector.cc
libRivetEvent_la_SOURCES = Event/Event.cc

libRivetTools_la_SOURCES = Logging.cc
libRivetTools_la_LDFLAGS = -L$(LOG4CPPLIBPATH) -R$(LOG4CPPLIBPATH)
libRivetTools_la_LIBADD  = libRivetMath.la libRivetEvent.la -l$(LOG4CPPLIBNAME)
