EXTRA_DIST = rivet.i ez_setup.py

rivet_wrap.cc rivet.py: rivet.i
	if test ! -e rivet_wrap.cc -o ! -e rivet.py; then \
      rm -f rivet_wrap.cc rivet.py; \
      swig -c++ -python -I$(top_srcdir)/include -o rivet_wrap.cc $<; \
      fi

if ENABLE_PYEXT

all-local: rivet_wrap.cc
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/*rivet*
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/lighthisto.py*

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
