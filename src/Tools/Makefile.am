noinst_LTLIBRARIES = libRivetMath.la libRivetEvent.la
pkglib_LTLIBRARIES = libRivetTools.la

libRivetMath_la_SOURCES  = Math/MathTools.cc Math/Vector.cc Math/Poincare.cc
libRivetEvent_la_SOURCES = Event/Event.cc

libRivetTools_la_SOURCES = Logging.cc
libRivetTools_la_LDFLAGS = -L$(LOG4CPPLIBPATH) -R$(LOG4CPPLIBPATH)
libRivetTools_la_LIBADD  = libRivetMath.la libRivetEvent.la -l$(LOG4CPPLIBNAME)
