SUBDIRS = Projections 
if ENABLE_ANALYSES
SUBDIRS += Analyses
endif
SUBDIRS += Tools . Test

lib_LTLIBRARIES  = libRivet.la

libRivet_la_SOURCES = \
  Event.cc Jet.cc \
  ProjectionApplier.cc Projection.cc \
  Analysis.cc AnalysisLoader.cc AnalysisInfo.cc \
  AnalysisHandler.cc Run.cc ProjectionHandler.cc HistoHandler.cc 

libRivet_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/src/Tools/yaml-cpp $(CPPFLAGS)

libRivet_la_LDFLAGS = -L$(prefix)/lib \
  -L$(GSLLIBPATH) -L$(GSLCBLASLIBPATH) -export-dynamic $(VERSIONINFOFLAGS)

FASTJETLIBADD = $(FASTJETCONFIGLIBADD)

libRivet_la_LIBADD  = \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la -ldl -lm \
  $(GSL_LDFLAGS) $(FASTJETLIBADD)
