noinst_LTLIBRARIES  = libRivetCore.la

libRivetCore_la_SOURCES = \
  Run.cc Event.cc Jet.cc Particle.cc ParticleName.cc \
  ProjectionApplier.cc Projection.cc \
  Analysis.cc AnalysisLoader.cc AnalysisInfo.cc \
  AnalysisHandler.cc ProjectionHandler.cc \
  Cuts.cc

if USE_BUILTIN_YAMLCPP
  SUBDIRS = yamlcpp
  libRivetCore_la_CPPFLAGS = -I$(srcdir)/yamlcpp/yaml-cpp/include -DYAMLCPP_API=3 -DYAML_NAMESPACE=RIVET_YAML $(AM_CPPFLAGS)
else
  libRivetCore_la_CPPFLAGS = $(AM_CPPFLAGS) $(YAMLCPP_CPPFLAGS) -DYAMLCPP_API=$(YAMLCPP_MAJOR_VERSION)
  libRivetCore_la_LDFLAGS = $(AM_LDFLAGS) $(YAMLCPP_LDFLAGS)
  libRivetCore_la_LIBADD = $(YAMLCPP_LIBS)
endif
