SUBDIRS = Projections Analyses Tools . Test


lib_LTLIBRARIES  = libRivet.la

libRivet_la_SOURCES = \
  Event.cc Jet.cc \
  ProjectionApplier.cc Projection.cc Analysis.cc AnalysisLoader.cc \
  AnalysisHandler.cc ProjectionHandler.cc HistoHandler.cc 

libRivet_la_CPPFLAGS = $(CPPFLAGS) $(AM_CPPFLAGS)

libRivet_la_LDFLAGS = \
  -L$(GSLLIBPATH) -L$(GSLCBLASLIBPATH) -export-dynamic $(VERSIONINFOFLAGS)

FASTJETLIBADD = $(FASTJETCONFIGLIBADD) \
  -lCDFConesPlugin -lJadePlugin -lD0RunIIConePlugin -lTrackJetPlugin \
  -lsiscone -lsiscone_spherical -lSISConePlugin

libRivet_la_LIBADD  = \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la -ldl -lm \
  -l$(GSLLIBNAME) -l$(GSLCBLASLIBNAME) \
  $(FASTJETLIBADD)
