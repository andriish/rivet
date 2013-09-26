SUBDIRS = rivet .

if ENABLE_PYEXT

PYEXT_ENV = CC="$(CC)" CXX="$(CXX)" ARCHFLAGS=""

## Always force setup.py, it's not good at guessing what needs to rebuild
all-local:
	$(PYEXT_ENV) $(PYTHON) setup.py install --install-lib=build/ --force

install-exec-local:
	$(PYEXT_ENV) $(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix) --force

uninstall-local:
	rm -rf $(DESTDIR)$(RIVET_PYTHONPATH)/rivet

clean-local:
	$(PYEXT_ENV) $(PYTHON) setup.py clean --all
	@rm -f $(builddir)/*.pyc
	@rm -rf $(builddir)/build
	@rm -rf dist

distclean-local:
	@rm -rf rivet.egg-info

endif
