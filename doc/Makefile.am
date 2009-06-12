## TODO: make variable names POSIX-compliant in rivet-manual.pdf target!

dist_noinst_SCRIPTS = mk-analysis-html mk-analysis-latex
EXTRA_DIST = compare-histos.txt heppennames.sty cone.png h-physrev3.bst hepnames.sty JHEP3.cls \
  make-plots.html bend.png hepnicenames.sty JHEP.bst make-plots.txt rivet-manual.tex preamble.tex \
  compare-histos.html hepparticles.sty maybemath.sty refs.bib rivet-manual.pdf thinker.png \
  hepunits.sty underscore.sty microtype.sty

DOCS = analyses.html

analyses.html: $(top_srcdir)/include/Rivet/Analyses mk-analysis-html
	LD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(LD_LIBRARY_PATH):$(prefix)/lib \
DYLD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(DYLD_LIBRARY_PATH):$(prefix)/lib \
./mk-analysis-html analyses.html


if ENABLE_PDFMANUAL
if WITH_PDFLATEX

DOCS += rivet-manual.pdf

analyses.tex: $(top_srcdir)/include/Rivet/Analyses mk-analysis-latex
	LD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(LD_LIBRARY_PATH):$(prefix)/lib \
DYLD_LIBRARY_PATH=$(top_srcdir)/src/.libs:$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(DYLD_LIBRARY_PATH):$(prefix)/lib \
./mk-analysis-latex analyses.tex

#LATEX	= pdflatex
LATEX	= pdflatex --interaction=nonstopmode
BIBTEX	= bibtex
MAKEINDEX = makeindex
RERUN = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined"
MAKEIDX = "^[^%]*\\makeindex"
COPY = if test -r $(<:%.tex=%.toc); then cp $(<:%.tex=%.toc) $(<:%.tex=%.toc.bak); fi
RM = rm -f

rivet-manual.pdf : rivet-manual.tex preamble.tex analyses.tex
	$(COPY);$(LATEX) $< 
	egrep $(MAKEIDX) $< && ($(MAKEINDEX) $(<:%.tex=%); $(COPY); $(LATEX) $<) > /dev/null; true
	egrep -c $(RERUNBIB) $(<:%.tex=%.log) && ($(BIBTEX) $(<:%.tex=%); $(COPY); $(LATEX) $<); true
	egrep $(RERUN) $(<:%.tex=%.log) && ($(COPY);$(LATEX) $<) > /dev/null; true
	egrep $(RERUN) $(<:%.tex=%.log) && ($(COPY);$(LATEX) $<) > /dev/null; true
	if cmp -s $(<:%.tex=%.toc) $(<:%.tex=%.toc.bak); then true; else $(LATEX) $< ; fi
	$(RM) $(<:%.tex=%.toc.bak)
	egrep -i "(Reference|Citation).*undefined" $(<:%.tex=%.log) ; true

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

all:
	@echo "Default make rule does nothing: use 'make doc'"

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
