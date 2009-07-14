EXTRA_DIST = rivet.i ez_setup.py

## Bundling SWIG-built sources with the tarball for SL4 compat
EXTRA_DIST += rivet_wrap.cc-HepMC23 rivet.py-HepMC23 \
              rivet_wrap.cc-HepMC24 rivet.py-HepMC24 \
              rivet_wrap.cc-HepMC25 rivet.py-HepMC25
rivet_wrap.cc rivet.py: rivet.i
	if test -e rivet_wrap.cc-HepMC$(HEPMC_VERSION); then cp rivet_wrap.cc-HepMC$(HEPMC_VERSION) rivet_wrap.cc; fi
	if test -e rivet.py-HepMC$(HEPMC_VERSION); then cp rivet.py-HepMC$(HEPMC_VERSION) rivet.py; fi
	if test ! -e rivet_wrap.cc -o ! -e rivet.py; then \
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
	@rm -rf $(builddir)/build
	@rm -rf dist
	@rm -f rivet.py
	@rm -f rivet_wrap.cc

distclean-local:
	@rm -rf rivet.egg-info

endif
