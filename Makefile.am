ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src pyext data include bin analyses test doc

EXTRA_DIST = GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-yoda=$(YODAPATH) --with-hepmc=$(HEPMCPATH) --with-fastjet=$(FASTJETPATH)

pkgconfigdir = $(libdir)/pkgconfig
dist_pkgconfig_DATA = rivet.pc

doc:
	cd doc && $(MAKE) doc

.PHONY : doc dox pyclean

clean-local:
	@rm -rf a.out


## Upload tarballs to HepForge
RSH = rsync
DEST = login.hepforge.org:rivet/downloads/
upload: dist
	$(RSH) $(DIST_ARCHIVES) $(DEST)
