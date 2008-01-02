if ENABLE_PYRIVET

SWIG_SOURCES  = pyrivet.i
BUILT_SOURCES = pyrivet_wrap.cc
pyrivet_wrap.cc : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) $(AM_CPPFLAGS) -o $@ $<

python_PYTHON = pyrivet.py
pyexec_LTLIBRARIES = _pyrivet.la
_pyrivet_la_SOURCES = SWIGTest.cc pyrivet_wrap.cc $(SWIG_SOURCES)
_pyrivet_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(HEPMCINCPATH) $(SWIG_PYTHON_CPPFLAGS) 
_pyrivet_la_CXXFLAGS = -Wno-long-long
_pyrivet_la_LDFLAGS = $(AM_LDFLAGS) -L$(HEPMCLIBPATH) -module
_pyrivet_la_LIBADD = -lHepMC -lRivet

all-local:
	@sh fixpyrivet

clean-local:
	@rm -f $(BUILT_SOURCES) pyrivet.py pyrivet.pyc

endif
