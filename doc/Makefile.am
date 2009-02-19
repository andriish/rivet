DOCNAME = rivet-manual

## TODO: make variable names POSIX-compliant!

LATEX	= pdflatex
BIBTEX	= bibtex
MAKEINDEX = makeindex
RERUN = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined"
MAKEIDX = "^[^%]*\\makeindex"

COPY = if test -r $(<:%.tex=%.toc); then cp $(<:%.tex=%.toc) $(<:%.tex=%.toc.bak); fi
RM = rm -f

.PHONY = doc

doc: $(DOCNAME).pdf compare-histos.html make-plots.html
	@true

analyses.tex: $(top_srcdir)/include/Rivet/Analyses mk-analysis-latex
	./mk-analysis-latex analyses.tex

$(DOCNAME).pdf : $(DOCNAME).tex analyses.tex
	$(COPY);$(LATEX) $<
	egrep $(MAKEIDX) $< && ($(MAKEINDEX) $(<:%.tex=%); $(COPY); $(LATEX) $<) > /dev/null; true
	egrep -c $(RERUNBIB) $(<:%.tex=%.log) && ($(BIBTEX) $(<:%.tex=%); $(COPY); $(LATEX) $<); true
	egrep $(RERUN) $(<:%.tex=%.log) && ($(COPY);$(LATEX) $<) > /dev/null; true
	egrep $(RERUN) $(<:%.tex=%.log) && ($(COPY);$(LATEX) $<) > /dev/null; true
	if cmp -s $(<:%.tex=%.toc) $(<:%.tex=%.toc.bak); then true; else $(LATEX) $< ; fi
	$(RM) $(<:%.tex=%.toc.bak)
	egrep -i "(Reference|Citation).*undefined" $(<:%.tex=%.log) ; true

mostlyclean-local:
	rm -rf *.aux *.log *.toc

clean-local:
	rm -rf $(DOCNAME).pdf

compare-histos.html: compare-histos.txt
	asciidoc -a toc compare-histos.txt

make-plots.html: make-plots.txt
	asciidoc -a toc make-plots.txt
