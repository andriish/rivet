SUBDIRS = Core Tools Projections AnalysisTools

lib_LTLIBRARIES = libRivet.la

libRivet_la_SOURCES =

libRivet_la_LDFLAGS = -export-dynamic -avoid-version -L$(YODALIBPATH)

libRivet_la_LIBADD = \
Core/libRivetCore.la \
Projections/libRivetProjections.la \
Tools/libRivetTools.la \
AnalysisTools/libRivetAnalysisTools.la \
-lYODA -ldl -lm  \
$(FASTJETLIBADD) -llwtnn

if ENABLE_HEPMC_3
libRivet_la_LIBADD += -lHepMC3 -lHepMC3search -llwtnn
libRivet_la_LDFLAGS += -L$(HEPMC3LIBPATH)
else
libRivet_la_LIBADD += -lHepMC -llwtnn
libRivet_la_LDFLAGS += -L$(HEPMCLIBPATH)
endif
