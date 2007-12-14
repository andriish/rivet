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
  $(FASTJETLIB) \
  -l$(GSLLIBNAME) -l$(GSLCBLASLIBNAME) \
  -l$(HEPPIDLIBNAME) -l$(HEPPDTLIBNAME) \
  -lSISConePlugin -lsiscone -lCDFConesPlugin
#  -lCDFcode -lCDFConesPlugin \
#  -lsiscone -lSISConePlugin
#  -l$(FASTJETLIBPATH)/libCDFcode -l$(FASTJETLIBPATH)/libCDFConesPlugin
#  -L$(FASTJETLIBPATH) #-lCDFcode -lCDFConesPlugin
if WITH_ROOT
  # libRivet_la_LIBADD += $(ROOTLIBS)	
  libRivet_la_LIBADD += $(ROOTLIBS)	
endif


## "rivetgun" executable
if ENABLE_RIVETGUN
  bin_PROGRAMS      = rivetgun
  rivetgun_SOURCES  = RivetgunInterface.cc RivetgunRunner.cc RivetgunCommandline.cc
  rivetgun_CPPFLAGS = $(AM_CPPFLAGS) -I$(AGILEINCPATH)
  rivetgun_LDFLAGS  = $(AM_LDFLAGS) -export-dynamic
  rivetgun_LDADD    = $(FCLIBS) $(AGILELIB) $(HEPMCLIB) $(GSLLIB) $(GSLCBLASLIB) libRivet.la
if WITH_ROOT
  rivetgun_CPPFLAGS += -DHAVE_ROOT $(ROOTCPPFLAGS)
endif
endif
