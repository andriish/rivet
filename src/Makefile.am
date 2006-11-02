SUBDIRS = Projections Analysis Tools Test

## Workaround for OS X, automake doesn't realize it's meant to use C++
#LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
#       $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

pkglib_LTLIBRARIES  = libRivet.la
libRivet_la_LDFLAGS = -version-info 1:0:0 -export-dynamic -L$(CLHEPLIBPATH)
libRivet_la_SOURCES = RivetHandler.cc RivetInfo.cc
libRivet_la_LIBADD  = Analysis/libRivetAnalysis.la \
                      Projections/libRivetProjections.la \
                      Tools/libRivetMath.la \
                      -l$(CLHEPLIBNAME) -ldl 
