EXTRA_DIST = fixpyrivet SWIGTest.cc rivet.py

if ENABLE_PYRIVET

# SWIG_SOURCES  = pyrivet.i
# BUILT_SOURCES = pyrivet_wrap.cc
# python_PYTHON = pyrivet.py

# pyrivet.py pyrivet_wrap.cc : $(SWIG_SOURCES)
# 	$(SWIG) $(SWIG_PYTHON_OPT) $(AM_CPPFLAGS) -o pyrivet_wrap.cc $<

# pyexec_LTLIBRARIES = _pyrivet.la
# _pyrivet_la_SOURCES = SWIGTest.cc pyrivet_wrap.cc $(SWIG_SOURCES)
# _pyrivet_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(HEPMCINCPATH) $(SWIG_PYTHON_CPPFLAGS) 
# _pyrivet_la_CXXFLAGS = -Wno-long-long
# _pyrivet_la_LDFLAGS = $(AM_LDFLAGS) -L$(HEPMCLIBPATH) -module
# _pyrivet_la_LIBADD = -lHepMC $(top_builddir)/src/libRivet.la

# all-local:
# 	@sh fixpyrivet

# clean-local:
# 	@rm -f $(BUILT_SOURCES) pyrivet.py pyrivet.pyc


# ## If not building Python extension
else 

EXTRA_DIST += pyrivet.i

endif
