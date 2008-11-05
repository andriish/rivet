bin_SCRIPTS = rivet-config 
dist_bin_SCRIPTS = aida2flat aida2root
EXTRA_DIST = flat2aida

if ENABLE_PYRIVET
dist_bin_SCRIPTS += rivet
#TESTS_ENVIRONMENT = RIVET_ANALYSIS_PATH=$(top_srcdir)/plugindemo
#TESTS = "rivet --list-analyses"
else
EXTRA_DIST += rivet
endif
