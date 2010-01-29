SUBDIRS = Core Projections Tools .
if ENABLE_ANALYSES
SUBDIRS += Analyses
endif

lib_LTLIBRARIES = libRivet.la

libRivet_la_SOURCES = 

libRivet_la_LDFLAGS = \
  -L$(prefix)/lib \
  -L$(GSLLIBPATH) -L$(GSLCBLASLIBPATH) \
  -export-dynamic $(VERSIONINFOFLAGS)

libRivet_la_LIBADD  = \
  Core/libRivetCore.la \
  Projections/libRivetProjections.la \
  Tools/libRivetTools.la \
  -ldl -lm \
  $(GSL_LDFLAGS) $(FASTJETCONFIGLIBADD)
