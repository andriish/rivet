SUBDIRS = Projections Analyses Tools . Test

## "libRivet" library
lib_LTLIBRARIES  = libRivet.la
libRivet_la_LDFLAGS = \
  -L$(HEPPIDLIBPATH) -L$(HEPPDTLIBPATH) \
  -L$(GSLLIBPATH) -L$(GSLCBLASLIBPATH) \
  -export-dynamic -version-info 1:0:0
libRivet_la_SOURCES = \
  Cuts.cc \
  Projection.cc \
  Analysis.cc \
  AnalysisLoader.cc \
  AnalysisHandler.cc 
libRivet_la_LIBADD  = \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la -ldl -lm \
  $(CLHEPLIB) $(KTJETLIB) $(FASTJETLIB) \
  -l$(GSLLIBNAME) -l$(GSLCBLASLIBNAME) \
  -l$(HEPPIDLIBNAME) -l$(HEPPDTLIBNAME)


## "rivet" executable
bin_PROGRAMS   = rivet
rivet_SOURCES  = RivetInterface.cc
rivet_CPPFLAGS = $(AM_CPPFLAGS)
rivet_LDFLAGS  = -export-dynamic
rivet_LDADD    = libRivet.la $(HEPMCLIB) $(CLHEPLIB) $(GSLLIB) $(GSLCBLASLIB)


## "rivetgun" executable
if ENABLE_RIVETGUN
  bin_PROGRAMS  += rivetgun
  rivetgun_SOURCES  = RivetGunInterfaceDynamic.cc RivetGunRunner.cc Tools/CommandlineRivetGun.cc
  rivetgun_CPPFLAGS = $(AM_CPPFLAGS) -I$(AGILEINCPATH)
  rivetgun_LDFLAGS  = $(AM_LDFLAGS) -export-dynamic
  rivetgun_LDADD    = $(FCLIBS) $(RIVETLIB) $(AGILELIB) $(HEPMCLIB)
endif


## Optional ROOT histogramming
if WITH_ROOT
  rivet_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS)
  libRivet_la_LIBADD += $(ROOTLIBS)	
endif
