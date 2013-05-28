SUBDIRS = rivet .

if ENABLE_PYEXT

PYEXT_ENV = CC=$(CC) CXX=$(CXX) ARCHFLAGS=""

all-local: rivet/rivetwrap_wrap.cc
	$(PYEXT_ENV) $(PYTHON) setup.py build

install-exec-local:
	$(PYEXT_ENV) $(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

## TODO: Really want DESTDIR here?
uninstall-local:
	rm -rf $(DESTDIR)$(RIVET_PYTHONPATH)/rivet
	rm -rf $(DESTDIR)$(RIVET_PYTHONPATH)/Rivet-*.egg-info

clean-local:
	$(PYEXT_ENV) $(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc
	@rm -rf $(builddir)/build
	@rm -rf dist

distclean-local:
	@rm -rf rivet.egg-info

endif
