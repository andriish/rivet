EXTRA_DIST = rivet.i

#if PYTHON_USE
if ENABLE_PYRIVET

all-local: rivet.i
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

uninstall-local:
	rm -rf $(libdir)/python*/site-packages/rivet.*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc

distclean-local:
	@rm -rf *.egg-info
	@rm -rf dist
	@rm -f rivet_wrap.* rivet.py*

endif
