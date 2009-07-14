EXTRA_DIST = rivet.i ez_setup.py

## Bundling SWIG-built sources with the tarball for SL4 compat
EXTRA_DIST += rivet_wrap.cc-HepMC24 rivet.py-HepMC24 \
              rivet_wrap.cc-HepMC25 rivet.py-HepMC25
rivet_wrap.cc rivet.py: rivet.i
	@cp rivet_wrap.cc-HepMC$(HEPMC_VERSION) rivet_wrap.cc
	@cp rivet.py-HepMC$(HEPMC_VERSION) rivet.py
	if test ! -e rivet_wrap.cc; then \
      rm -f rivet_wrap.cc rivet.py; \
      swig -c++ -python -I$(HEPMCINCPATH) -I$(top_srcdir)/include -o rivet_wrap.cc $<; \
      fi

if ENABLE_PYRIVET

all-local: rivet_wrap.cc
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/rivet*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc
	@rm -f $(builddir)/build
	@rm -rf dist

distclean-local:
	@rm -rf rivet.egg-info
	@rm -f rivet.py
	@rm -f rivet_wrap.cc

endif
