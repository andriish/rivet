lib_LTLIBRARIES = RivetExampleAnalyses.la

RivetExampleAnalyses_la_SOURCES = \
    ExampleAnalysis.cc \
    ExampleTree.cc

# ## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -L$(FASTJETLIBPATH)
LIBS = -lCDFConesPlugin -lSISConePlugin -lsiscone -lJadePlugin
