ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src pyext data include bin plugindemo

doc:
	cd doc && $(MAKE) doc

.PHONY : doc dox

clean-local:
	@rm -rf a.out

## Doxygen
if WITH_DOXYGEN
dox:
	$(DOXYGEN) $(top_srcdir)/Doxyfile || true
install-dox:
	$(install_sh_DATA) doxy/html/* $(DESTDIR)$(pkgdatadir)/doc/html
#all-local: dox
#	@true
#install-data-local:
#	$(install_sh_DATA) doxy/html/* $(DESTDIR)$(pkgdatadir)/doc/html
mostlyclean-local:
	@rm -rf $(top_builddir)/doxy
uninstall-local:
	@rm -rf $(DESTDIR)$(pkgdatadir)/doc; \
	test ! -d $(DESTDIR)$(pkgdatadir) || rmdir --ignore-fail-on-non-empty $(DESTDIR)$(pkgdatadir)
endif

## Remove SVN dirs
dist-hook:
	@rm -rf `find $(distdir) -name ".svn"`
