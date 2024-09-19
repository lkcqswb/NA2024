OBJ = answer.tex

doc: ${OBJ}
	@pdflatex ${OBJ}
	@bibtex main
	@pdflatex ${OBJ}
	@pdflatex ${OBJ}

.PHONY: clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz *.toc *.bbl *.run.xml *.blg