bin_SCRIPTS = rivet-config 

dist_bin_SCRIPTS = \
	aida2flat aida2root flat2aida \
	compare-histos make-plots
EXTRA_DIST = 

RIVETPROGS = \
    rivet \
    rivet-mkanalysis rivet-buildplugin \
    rivet-chopbins rivet-rmgaps rivet-rescale \
    rivet-mergeruns rivet-mkhtml

if ENABLE_PYEXT
dist_bin_SCRIPTS += $(RIVETPROGS)
else
EXTRA_DIST += $(RIVETPROGS)
endif
