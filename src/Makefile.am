SUBDIRS = Core Tools Projections
if ENABLE_ANALYSES
SUBDIRS += Analyses
ANA_LIBADD = Analyses/libRivetAnalysisTools.la
endif
SUBDIRS += .

lib_LTLIBRARIES = libRivet.la

libRivet_la_SOURCES =

libRivet_la_LDFLAGS = -export-dynamic -avoid-version \
  Core/libRivetCore.la Projections/libRivetProjections.la Tools/libRivetTools.la $(ANA_LIBADD) \
  -L$(YODALIBPATH) -lYODA -L$(HEPMCLIBPATH) -lHepMC -ldl -lm  $(GSL_LDFLAGS) $(FASTJETCONFIGLIBADD)

if USE_BUILTIN_YAMLCPP
libRivet_la_LDFLAGS += -L$(srcdir)/Core/yamlcpp/yaml-cpp -lrivet-yaml-cpp
endif
