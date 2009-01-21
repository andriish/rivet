## TODO: make the doc rule run until a successful build (or > 3 runs)
doc: rivet-manual.pdf
	@true

rivet-manual.pdf: rivet-manual.tex
	while true; do pdflatex rivet-manual.tex; \
    if grep -i rerun rivet-manual.log; then true; else break; fi; done 

mostlyclean-local:
	rm -rf *.aux *.log *.toc

clean-local:
	rm -rf rivet-manual.pdf
