## TODO: make variable names POSIX-compliant in rivet-manual.pdf target!

dist_noinst_SCRIPTS = mk-analysis-html mk-analysis-latex
EXTRA_DIST = compare-histos.txt heppennames.sty cone.png h-physrev3.bst hepnames.sty JHEP3.cls \
  make-plots.html bend.png hepnicenames.sty JHEP.bst make-plots.txt rivet-manual.tex preamble.tex \
  compare-histos.html hepparticles.sty maybemath.sty refs.bib rivet-manual.pdf thinker.png \
  hepunits.sty underscore.sty microtype.sty

DOCS = analyses.html

if ENABLE_PYEXT

analyses.html: $(top_srcdir)/include/Rivet/Analyses mk-analysis-html
	LD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(LD_LIBRARY_PATH):$(prefix)/lib \
DYLD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(DYLD_LIBRARY_PATH):$(prefix)/lib \
RIVET_DATA_PATH=$(top_srcdir)/data/anainfo \
./mk-analysis-html analyses.html

endif


if ENABLE_PDFMANUAL
if WITH_PDFLATEX

DOCS += rivet-manual.pdf


if ENABLE_PYEXT

analyses.tex: $(top_srcdir)/include/Rivet/Analyses mk-analysis-latex
	LD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(LD_LIBRARY_PATH):$(prefix)/lib \
DYLD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(DYLD_LIBRARY_PATH):$(prefix)/lib \
RIVET_DATA_PATH=$(top_srcdir)/data/anainfo \
./mk-analysis-latex analyses.tex

else

analyses.tex:
	> analyses.tex

endif


#LATEX	= pdflatex
LATEX	= pdflatex --interaction=nonstopmode
BIBTEX	= bibtex
MAKEINDEX = makeindex
RERUN = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined"
MAKEIDX = "^[^%]*\\makeindex"
RM = rm -f

DOCNAME = rivet-manual
$(DOCNAME).pdf : $(DOCNAME).tex preamble.tex analyses.tex
	$(LATEX) $<; true
	egrep $(MAKEIDX) $< && ($(MAKEINDEX) $(DOCNAME) && cp $(DOCNAME).toc $(DOCNAME).toc.bak && $(LATEX) $<) > /dev/null; true
	egrep -c $(RERUNBIB) $(DOCNAME).log && ($(BIBTEX) $(DOCNAME) && cp $(DOCNAME).toc $(DOCNAME).toc.bak && $(LATEX) $<); true
	egrep $(RERUN) $(DOCNAME).log && cp $(DOCNAME).toc $(DOCNAME).toc.bak && $(LATEX) $<; true
	egrep $(RERUN) $(DOCNAME).log && cp $(DOCNAME).toc $(DOCNAME).toc.bak && $(LATEX) $<; true
	if cmp -s $(DOCNAME).toc $(DOCNAME).toc.bak; then true; else $(LATEX) $<; true; fi
	$(RM) $(DOCNAME).toc.bak; true

endif
endif


if WITH_ASCIIDOC

DOCS += compare-histos.html make-plots.html 

compare-histos.html: compare-histos.txt
	asciidoc -a toc compare-histos.txt

make-plots.html: make-plots.txt
	asciidoc -a toc make-plots.txt

endif


################


.PHONY = all doc upload

#all:
#	@echo "Default make rule does nothing: use 'make doc'"

doc: $(DOCS)
	@true

## TODO: Put Rivet version string in PDF filename for upload?
RSH=rsync
DEST=login.hepforge.org:rivet/public_html/
upload: $(DOCS)
	$(RSH) $? $(DEST)

mostlyclean-local:
	rm -rf *.aux *.log *.toc

clean-local:
	rm -rf $(DOCS)


## Install!
pkgdata_DATA = $(DOCS)
