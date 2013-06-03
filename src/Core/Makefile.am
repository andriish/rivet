noinst_LTLIBRARIES  = libRivetCore.la

libRivetCore_la_SOURCES = \
  Event.cc Jet.cc Particle.cc ParticleName.cc \
  ProjectionApplier.cc Projection.cc \
  Analysis.cc AnalysisLoader.cc AnalysisInfo.cc \
  AnalysisHandler.cc Run.cc ProjectionHandler.cc \
  Cuts.cc

libRivetCore_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(YAML_CPPINCPATH) -DYAMLCPP_API_VERSION=$(YAMLCPP_MAJOR_VERSION)
