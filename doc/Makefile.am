## TODO: make variable names POSIX-compliant in rivet-manual.pdf target!


DOCS = analyses.html

analyses.html: $(top_srcdir)/include/Rivet/Analyses mk-analysis-html
	./mk-analysis-html analyses.html


if WITH_PDFLATEX

DOCS += rivet-manual.pdf

analyses.tex: $(top_srcdir)/include/Rivet/Analyses mk-analysis-latex
	./mk-analysis-latex analyses.tex

LATEX	= pdflatex
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

RSH=rsync
DEST=login.hepforge.org:rivet/public_html/
upload: $(DOCS)
	$(RSH) $? $(DEST)

mostlyclean-local:
	rm -rf *.aux *.log *.toc

clean-local:
	rm -rf $(DOCS)


## Install!
#pkgdata_DATA = $(DOCS)
