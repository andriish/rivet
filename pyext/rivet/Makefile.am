EXTRA_DIST = \
	__init__.py \
    hepdatautils.py \
    spiresbib.py util.py \
    plotinfo.py aopaths.py \
    core.pyx rivet.pxd core.cpp \
	hepdatapatches


MAINTAINERCLEANFILES = core.cpp
BUILT_SOURCES = core.cpp

if WITH_CYTHON
.NOTPARALLEL:

RIVETINCLUDE = $(top_srcdir)/include/Rivet/
core.cpp: core.pyx rivet.pxd $(RIVETINCLUDE)/Analysis.hh $(RIVETINCLUDE)/AnalysisHandler.hh $(RIVETINCLUDE)/AnalysisLoader.hh $(RIVETINCLUDE)/Run.hh
	$(CYTHON) $(srcdir)/core.pyx --cplus -2 \
	   -I $(srcdir) -I $(srcdir)/include \
	   -I $(builddir) -I $(builddir)/include \
	   -o $@
else
core.cpp:
	@echo "Not (re)generating core.cpp since Cython is not installed"
endif
