SUBDIRS = Projections Analysis Tools 
#Test

## Workaround for OS X, automake doesn't realize it's meant to use C++
#LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
#       $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

pkglib_LTLIBRARIES  = libRivet.la
libRivet_la_LDFLAGS = -version-info 1:0:0 -export-dynamic \
                      -L$(CLHEPLIBPATH)   -R$(CLHEPLIBPATH) \
                      -L$(LOG4CPPLIBPATH) -R$(LOG4CPPLIBPATH) \
                      -L$(HEPPIDLIBPATH)  -R$(HEPPIDLIBPATH)
libRivet_la_SOURCES = RivetHandler.cc RivetInfo.cc Logging.cc
libRivet_la_LIBADD  = Analysis/libRivetAnalysis.la \
                      Projections/libRivetProjections.la \
                      Tools/libRivetMath.la -ldl \
                      -l$(CLHEPLIBNAME) -l$(LOG4CPPLIBNAME) -l$(HEPPIDLIBNAME)

SUBDIRS += Test
