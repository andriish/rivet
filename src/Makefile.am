SUBDIRS = Projections Analysis Tools . Test

## Workaround for OS X, automake doesn't realize it's meant to use C++
#LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) \
#       $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

lib_LTLIBRARIES  = libRivet.la
libRivet_la_LDFLAGS = -version-info 1:0:0 -export-dynamic \
                      -L$(CLHEP_VECTORLIBPATH) -R$(CLHEP_VECTORLIBPATH) \
                      -L$(CLHEP_MATRIXLIBPATH) -R$(CLHEP_MATRIXLIBPATH) \
                      -L$(HEPPIDLIBPATH)  -R$(HEPPIDLIBPATH) \
                      -L$(HEPPDTLIBPATH)  -R$(HEPPDTLIBPATH) \
                      -L$(KTJETLIBPATH) -R$(KTJETLIBPATH)
libRivet_la_SOURCES = RivetHandler.cc RivetInfo.cc
libRivet_la_LIBADD  = Analysis/libRivetAnalysis.la \
                      Projections/libRivetProjections.la \
                      Tools/libRivetTools.la -ldl \
					  -l$(CLHEP_MATRIXLIBNAME) -l$(CLHEP_VECTORLIBNAME) \
                      -l$(HEPPIDLIBNAME) -l$(HEPPDTLIBNAME) \
                      -l$(KTJETLIBNAME)

bin_PROGRAMS = rivet
rivet_SOURCES = RivetInterface.cc
rivet_LDFLAGS = -L$(HEPMCLIBPATH) -R$(HEPMCLIBPATH)
rivet_LDADD = -lRivet -l$(HEPMCLIBNAME)
