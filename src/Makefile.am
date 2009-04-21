SUBDIRS = Projections Analyses Tools . Test


lib_LTLIBRARIES  = libRivet.la

libRivet_la_SOURCES = \
  Jet.cc \
  ProjectionApplier.cc Projection.cc Analysis.cc AnalysisLoader.cc \
  AnalysisHandler.cc ProjectionHandler.cc HistoHandler.cc 

libRivet_la_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS)

libRivet_la_LDFLAGS = \
  -L$(GSLLIBPATH) -L$(GSLCBLASLIBPATH) -export-dynamic $(VERSIONINFOFLAGS)

FASTJETLIBADD = $(FASTJETCONFIGLIBADD) \
  -lCDFConesPlugin -lsiscone -lSISConePlugin
if ENABLE_JADE
FASTJETLIBADD += -lJadePlugin
endif

libRivet_la_LIBADD  = \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la -ldl -lm \
  -l$(GSLLIBNAME) -l$(GSLCBLASLIBNAME) \
  $(FASTJETLIBADD)
