SUBDIRS = Projections Analyses Tools . Test

lib_LTLIBRARIES  = libRivet.la
libRivet_la_LDFLAGS = \
  -L$(HEPPIDLIBPATH) -L$(HEPPDTLIBPATH) \
  -export-dynamic -version-info 1:0:0
libRivet_la_SOURCES = \
  Particle.cc \
  Cuts.cc \
  Projection.cc \
  Analysis.cc \
  AnalysisLoader.cc \
  AnalysisHandler.cc 
libRivet_la_LIBADD  = \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la -ldl \
  $(CLHEPLIB) $(KTJETLIB) \
  -l$(HEPPIDLIBNAME) -l$(HEPPDTLIBNAME) 


bin_PROGRAMS = rivet
rivet_SOURCES = RivetInterface.cc
rivet_LDFLAGS = -L$(top_srcdir)/src -export-dynamic
rivet_LDADD = -lRivet $(HEPMCLIB) $(CLHEPLIB) 
rivet_CPPFLAGS = $(AM_CPPFLAGS)


if WITH_ROOT
  rivet_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS)
  libRivet_la_LIBADD += $(ROOTLIBS)	
endif

if WITH_FASTJET
  rivet_CPPFLAGS += -DHAVE_FASTJET
  ## TODO: Use abs path rather than -l flag
  libRivet_la_LIBADD += -l$(FASTJETLIBNAME)
endif

