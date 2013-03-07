bin_SCRIPTS = rivet-config
dist_bin_SCRIPTS = make-plots
EXTRA_DIST =

RIVETPROGS = \
    rivet \
    rivet-mkanalysis rivet-buildplugin \
	rivet-findid rivet-which \
    rivet-cmphistos rivet-mkhtml
    # TODO: reinstate  rivet-rescale rivet-mergeruns rivet-merge-CDF_2012_NOTE10874
if ENABLE_PYEXT
dist_bin_SCRIPTS += $(RIVETPROGS)
else
EXTRA_DIST += $(RIVETPROGS)
endif

noinst_PROGRAMS = rivet-nopy
rivet_nopy_SOURCES = rivet-nopy.cc
rivet_nopy_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)
rivet_nopy_LDADD = $(top_builddir)/src/libRivet.la -lHepMC
rivet_nopy_LDFLAGS = -L$(HEPMCLIBPATH)
