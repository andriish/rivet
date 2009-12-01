noinst_LTLIBRARIES  = libRivetCore.la

libRivetCore_la_SOURCES = \
  Event.cc Jet.cc Particle.cc \
  ProjectionApplier.cc Projection.cc \
  Analysis.cc AnalysisLoader.cc AnalysisInfo.cc \
  AnalysisHandler.cc Run.cc ProjectionHandler.cc HistoHandler.cc

libRivetCore_la_CPPFLAGS = \
  $(AM_CPPFLAGS) -I$(top_srcdir)/src/Tools/yaml-cpp $(CPPFLAGS)
