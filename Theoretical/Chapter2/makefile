

report: doc

doc: answer.tex
	@pdflatex answer.tex
	@bibtex answer
	@pdflatex answer.tex
	@pdflatex answer.tex


.PHONY: clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz *.toc *.bbl *.run.xml *.blg *.bib