install-data-local:
	for i in `find $(srcdir) -name \*.mplstyle` ; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkgdatadir)/stylelib/$${i#$(srcdir)/}; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(pkgdatadir)/texmf

