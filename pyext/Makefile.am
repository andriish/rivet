EXTRA_DIST = rivet.i ez_setup.py lighthisto.py

rivet_wrap.cc rivet.py: rivet.i
	swig -c++ -python -I$(top_srcdir)/include -o rivet_wrap.cc $<;

if ENABLE_PYEXT

all-local: rivet_wrap.cc
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(prefix)/lib*/python*/site-packages/*rivet*
	rm -rf $(DESTDIR)$(prefix)/lib*/python*/site-packages/lighthisto.py*

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
