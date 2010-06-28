EXTRA_DIST = ez_setup.py lighthisto.py spiresbib.py

if ENABLE_PYEXT

SUBDIRS = rivet .

all-local: rivet/rivetwrap_wrap.cc
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(prefix)/lib*/python*/site-packages/rivet
	rm -rf $(DESTDIR)$(prefix)/lib*/python*/site-packages/lighthisto.py*
	rm -rf $(DESTDIR)$(prefix)/lib*/python*/site-packages/spiresbib.py*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc
	@rm -rf $(builddir)/build
	@rm -rf dist

distclean-local:
	@rm -rf rivet.egg-info

endif
