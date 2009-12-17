bin_SCRIPTS = rivet-config rivet-buildplugin
dist_bin_SCRIPTS = aida2flat aida2root compare-histos make-html make-plots rivet-mkanalysis rivet-chop-bins rivet-normalise-histos
EXTRA_DIST = flat2aida

if ENABLE_PYEXT
dist_bin_SCRIPTS += rivet
#TESTS_ENVIRONMENT = RIVET_ANALYSIS_PATH=$(top_srcdir)/plugindemo
#TESTS = "rivet --list-analyses"
else
EXTRA_DIST += rivet
endif
