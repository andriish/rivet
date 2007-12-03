ACLOCAL_AMFLAGS = -I acincludes

include_HEADERS = HepMCView.hh

bin_PROGRAMS = hepmclist
hepmclist_CPPFLAGS = -I$(HOME)/heplocal/include
hepmclist_SOURCES  = HepMCView.cc HepMCViewMain.cc
hepmclist_LDFLAGS = -L$(HOME)/heplocal/lib
hepmclist_LDADD = -lHepMC

BUILT_SOURCES = hepmcview_wrap.cc
SWIG_SOURCES  = hepmcview.i

pkgpython_PYTHON = hepmcview.py
pkgpyexec_LTLIBRARIES = _hepmcview.la
_hepmcview_la_SOURCES = HepMCView.cc hepmcview_wrap.cc $(SWIG_SOURCES)
_hepmcview_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(HOME)/heplocal/include
_hepmcview_la_LDFLAGS = -module -L$(HOME)/heplocal/lib
_hepmcview_la_LIBADD = -lHepMC

hepmcview_wrap.cc : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(HOME)/heplocal/include -o $@ $<
