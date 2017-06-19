AUTOMAKE_OPTIONS = -Wno-portability

## Flags for building all plugins
# TODO: remove?
AM_LDFLAGS = $(LDFLAGS) -module -avoid-version -L$(FASTJETLIBPATH)
LIBS = $(FASTJETCONFIGLIBADD)
lib_LTLIBRARIES =
noinst_LTLIBRARIES = libRivetAnalysisTools.la
libRivetAnalysisTools_la_SOURCES = \
    MC_ParticleAnalysis.cc \
    MC_JetAnalysis.cc \
    MC_JetSplittings.cc


if ENABLE_ANALYSES

# TODO: ??? remove?
BUILT_SOURCES =
DISTCLEANFILES =

PLUGIN_DIRS = $(wildcard plugin*)
PLUGIN_NAMES = $(subst plugin,,$(PLUGIN_DIRS))
PLUGIN_LIBS = $(patsubst plugin%,Rivet%Analyses.so,$(PLUGIN_DIRS))
PLUGIN_DATAFILES = $(wildcard plugin*/*.info plugin*/*.plot plugin*/*.yoda)

CLEANFILES = $(PLUGIN_LIBS)
EXTRA_DIST = $(PLUGIN_DIRS)

%.so:
	RIVET_BUILDPLUGIN_BEFORE_INSTALL=1 bash $(top_srcdir)/bin/rivet-buildplugin $@ $^

RivetALICEAnalyses.so: $(wildcard pluginALICE/*.cc)
RivetATLASAnalyses.so: $(wildcard pluginATLAS/*.cc)
RivetCDFAnalyses.so: $(wildcard pluginCDF/*.cc)
RivetCMSAnalyses.so: $(wildcard pluginCMS/*.cc)
RivetD0Analyses.so: $(wildcard pluginD0/*.cc)
RivetHERAAnalyses.so: $(wildcard pluginHERA/*.cc)
RivetLEPAnalyses.so: $(wildcard pluginLEP/*.cc)
RivetLHCbAnalyses.so: $(wildcard pluginLHCb/*.cc)
RivetLHCfAnalyses.so: $(wildcard pluginLHCf/*.cc)
RivetMCAnalyses.so: $(wildcard pluginMC/*.cc)
RivetMiscAnalyses.so: $(wildcard pluginMisc/*.cc)
RivetPetraAnalyses.so: $(wildcard pluginPetra/*.cc)
RivetRHICAnalyses.so: $(wildcard pluginRHIC/*.cc)
RivetSPSAnalyses.so: $(wildcard pluginSPS/*.cc)
RivetTOTEMAnalyses.so: $(wildcard pluginTOTEM/*.cc)

all-local: $(PLUGIN_LIBS)

# TODO: need dist and install hooks to strip out .hg dirs

install-exec-local: $(PLUGIN_LIBS)
	$(INSTALL) $^ $(DESTDIR)$(libdir)

install-data-local: $(PLUGIN_DATAFILES)
	@echo $^
	$(install_sh_DATA) $^ $(DESTDIR)$(libdir)

uninstall-local:
	@cd $(DESTDIR)$(libdir) && rm -rf $(PLUGIN_LIBS)

endif
