EXTRA_DIST = rivet.i ez_setup.py

## Bundling SWIG-built sources with the tarball for SL4 compat
EXTRA_DIST += rivet_wrap.cc rivet.py
rivet_wrap.cc rivet.py: rivet.i
	swig -c++ -python -I$(HEPMCINCPATH) -I$(top_srcdir)/include -o rivet_wrap.cc $<


#if PYTHON_USE
if ENABLE_PYRIVET

all-local: rivet_wrap.cc
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(libdir)/python*/site-packages/rivet*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc

distclean-local:
	@rm -rf rivet.egg-info
	@rm -rf dist
	@rm -f rivet.py
	@rm -f rivet_wrap.cc

endif
