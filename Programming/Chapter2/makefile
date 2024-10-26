CXX = g++
CXXFLAGS = -std=c++11 -Wall
TARGETS = ProblemB ProblemC ProblemD ProblemE ProblemF

.PHONY: all compile run report clean

all: compile

compile: $(TARGETS)

$(TARGETS): %: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run: compile
	@for target in $(TARGETS); do \
		./$$target; \
	done
	@echo "All targets have been executed."
	$(MAKE) clean 

run_test:
	@for target in $(TARGETS); do \
		echo "Running $$target..."; \
		./$$target > $$target.txt; \
	done

report: compile run_test report.pdf clean clean_picture

report.pdf: report.tex
	@xelatex report.tex
	@bibtex report
	@xelatex report.tex
	@xelatex report.tex

clean:
	rm -f $(TARGETS) report.aux report.log report.bbl report.blg report-blx.bib report.run.xml ProblemB ProblemC ProblemD ProblemE ProblemF ProblemB.txt ProblemC.txt ProblemD.txt ProblemE.txt ProblemF.txt report.synctex.gz

clean_picture: rm -f ProblemCplot.png ProblemBplot.png ProblemEplot.png ProblemF_m=10_plot.png ProblemF_m=40_plot.png ProblemF_m=60_plot.png