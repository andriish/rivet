SUBDIRS = Core Tools Projections AnalysisTools

##BinnedHisto if ENABLE_ANALYSES
##BinnedHisto ANA_LIBADD += Analyses/libRivetAnalysisTools.la
##BinnedHisto endif

lib_LTLIBRARIES = libRivet.la

libRivet_la_SOURCES =

libRivet_la_LDFLAGS = -export-dynamic -avoid-version -L$(YODALIBPATH) -L$(HEPMCLIBPATH)

libRivet_la_LIBADD = \
Core/libRivetCore.la \
Projections/libRivetProjections.la \
Tools/libRivetTools.la \
AnalysisTools/libRivetAnalysisTools.la \
-lYODA -lHepMC -ldl -lm  \
$(FASTJETCONFIGLIBADD) \
$(GSL_LDFLAGS)
