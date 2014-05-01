EXTRA_DIST = \
	__init__.py util.py \
    plotinfo.py spiresbib.py \
    core.pyx rivet.pxd \
    core.cpp


if WITH_CYTHON
RIVETINCLUDE = $(top_srcdir)/include/Rivet/
core.cpp: core.pyx rivet.pxd $(RIVETINCLUDE)/Analysis.hh $(RIVETINCLUDE)/AnalysisHandler.hh $(RIVETINCLUDE)/AnalysisLoader.hh $(RIVETINCLUDE)/Run.hh
	cython core.pyx --cplus
else
core.cpp:
	@echo "Not (re)generating core.cpp since Cython is not installed"
endif


all-local: core.cpp
	@true

clean-local:
	rm -f core.cpp
