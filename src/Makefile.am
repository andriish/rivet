SUBDIRS = Projections Analyses Tools . Test

## "libRivet" library
lib_LTLIBRARIES  = libRivet.la

libRivet_la_LDFLAGS = \
  -L$(GSLLIBPATH) -L$(GSLCBLASLIBPATH) -export-dynamic $(VERSIONINFOFLAGS)

libRivet_la_SOURCES = \
  Cuts.cc Jet.cc \
  ProjectionApplier.cc Projection.cc Analysis.cc AnalysisLoader.cc \
  AnalysisHandler.cc ProjectionHandler.cc HistoHandler.cc 

libRivet_la_LIBADD  = \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la -ldl -lm \
  $(FASTJETCONFIGLIBADD) \
  -lCDFConesPlugin -lsiscone -lSISConePlugin \
  -l$(GSLLIBNAME) -l$(GSLCBLASLIBNAME)

libRivet_la_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS)

if WITH_ROOT
  libRivet_la_LIBADD += $(ROOTLIBS)
  libRivet_la_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS)
endif

if ENABLE_JADE
  libRivet_la_LIBADD += -lJadePlugin
endif


## "rivetgun" executable
EXTRA_PROGRAMS = rivetgun

if ENABLE_RIVETGUN
  bin_PROGRAMS      = rivetgun
  rivetgun_SOURCES  = RivetgunInterface.cc RivetgunRunner.cc RivetgunCommandline.cc
  rivetgun_CPPFLAGS = $(AM_CPPFLAGS) -I$(AGILEINCPATH) -DRIVETVERSION=\"$(PACKAGE_VERSION)\"
  rivetgun_LDFLAGS  = $(AM_LDFLAGS) -export-dynamic
  rivetgun_LDADD    = $(FCLIBS) $(AGILELIB) $(HEPMCLIB) $(GSLLIB) $(GSLCBLASLIB) libRivet.la
if WITH_ROOT
  rivetgun_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS)
endif
endif
