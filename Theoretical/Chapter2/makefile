

report: doc

doc: answer.tex
	@pdflatex HW2_3220103383_LuoKaicheng.tex
	@bibtex HW2_3220103383_LuoKaicheng
	@pdflatex HW2_3220103383_LuoKaicheng.tex
	@pdflatex HW2_3220103383_LuoKaicheng.tex


.PHONY: clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz *.toc *.bbl *.run.xml *.blg *.bib