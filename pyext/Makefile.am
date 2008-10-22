EXTRA_DIST = fixpyrivet

if ENABLE_PYRIVET

SWIG_SOURCES  = rivet.i
BUILT_SOURCES = rivet_wrap.cc
python_PYTHON = rivet.py

rivet.py rivet_wrap.cc : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) $(AM_CPPFLAGS) -o rivet_wrap.cc $<

pyexec_LTLIBRARIES = _rivet.la
_rivet_la_SOURCES = rivet_wrap.cc $(SWIG_SOURCES)
_rivet_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(HEPMCINCPATH) $(SWIG_PYTHON_CPPFLAGS) 
_rivet_la_CXXFLAGS = -Wno-long-long
_rivet_la_LDFLAGS = $(AM_LDFLAGS) -L$(HEPMCLIBPATH) -module
_rivet_la_LIBADD = -lHepMC $(top_builddir)/src/libRivet.la

all-local:
	@sh fixpyrivet

clean-local:
	@rm -f $(BUILT_SOURCES) rivet.py rivet.pyc


## If not building Python extension
else 

EXTRA_DIST += rivet.i

endif
