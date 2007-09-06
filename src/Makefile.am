SUBDIRS = Projections Analysis Tools . Test

lib_LTLIBRARIES  = libRivet.la
libRivet_la_LDFLAGS = \
  -version-info 1:0:0 -export-dynamic \
  -L$(HEPPIDLIBPATH) -L$(HEPPDTLIBPATH) 
libRivet_la_SOURCES = AnalysisHandler.cc Particle.cc Cuts.cc
libRivet_la_LIBADD  = \
  Analysis/libRivetAnalysis.la \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la -ldl \
  $(CLHEPLIB) $(KTJETLIB) \
  -l$(HEPPIDLIBNAME) -l$(HEPPDTLIBNAME) 


bin_PROGRAMS = rivet
rivet_SOURCES = RivetInterface.cc
rivet_LDFLAGS = -L$(top_srcdir)/src
rivet_LDADD = -lRivet $(HEPMCLIB) 
rivet_CPPFLAGS = $(AM_CPPFLAGS) -I$(RIVETINCPATH)
if WITH_ROOT
  rivet_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS)
  libRivet_la_LIBADD += $(ROOTLIBS)	
endif
if WITH_FASTJET
  rivet_CPPFLAGS += -DHAVE_FASTJET
  libRivet_la_LIBADD += $(FASTJETLIB)	
endif

check-local:
	echo "HepPDT = $(HEPPDTLIB)"
