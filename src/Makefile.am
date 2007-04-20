SUBDIRS = Projections Analysis Tools . Test

lib_LTLIBRARIES  = libRivet.la
libRivet_la_LDFLAGS = -version-info 1:0:0 -export-dynamic \
                      -L$(CLHEPLIBPATH) -R$(CLHEPLIBPATH) \
                      -L$(HEPPIDLIBPATH)  -R$(HEPPIDLIBPATH) \
                      -L$(HEPPDTLIBPATH)  -R$(HEPPDTLIBPATH) \
                      -L$(KTJETLIBPATH) -R$(KTJETLIBPATH)
libRivet_la_SOURCES = AnalysisHandler.cc Particle.cc Cut.cc
libRivet_la_LIBADD  = Analysis/libRivetAnalysis.la \
                      Projections/libRivetProjections.la \
                      Tools/libRivetTools.la -ldl \
					  -l$(CLHEPLIBNAME)  \
                      -l$(HEPPIDLIBNAME) \
                      -l$(HEPPDTLIBNAME) \
                      -l$(KTJETLIBNAME)

bin_PROGRAMS = rivet
rivet_SOURCES = RivetInterface.cc
rivet_LDFLAGS = -L$(top_srcdir/src) -L$(HEPMCLIBPATH) -R$(HEPMCLIBPATH)
rivet_LDADD = -lRivet -l$(HEPMCLIBNAME)
