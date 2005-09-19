BUILT_SOURCES = done-all-links
CLEANFILES = done-all-links \
Rivet/Projections \
Rivet/Analysis \
Rivet/Config \
Rivet/CLHEPWrap

my_directories = Config CLHEPWrap Projections Analysis

done-all-links:
@EMPTY@ifdef SHOWCOMMAND
	mkdir -p Rivet
	cd Rivet; \
	for i in $(my_directories); \
	do \
	$(LN_S) -f ../$(top_srcdir)/$$i ; \
	done
	echo "stamp" > done-all-links
@EMPTY@else
	@echo "sym-linking header files..."
	@mkdir -p Rivet
	@cd Rivet; \
	for i in $(my_directories); \
	do \
	$(LN_S) -f ../$(top_srcdir)/$$i ; \
	done
	@echo "stamp" > done-all-links
@EMPTY@endif

install-data-local:
	find Rivet -follow -name '*.h' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;
	find Rivet -follow -name '*.icc' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;
	find Rivet -follow -name '*.tcc' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;
	find Rivet -follow -name '*.fh'  \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/Rivet
