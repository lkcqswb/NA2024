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

run_test:
	@for target in $(TARGETS); do \
		echo "Running $$target..."; \
		./$$target > $$target.txt; \
	done

report: compile run_test report.pdf clean

report.pdf: report.tex
	@pdflatex report.tex
	@bibtex report
	@pdflatex report.tex
	@pdflatex report.tex

clean:
	rm -f $(TARGETS) report.aux report.log report.bbl report.blg report-blx.bib report.run.xml ProblemB ProblemC ProblemD ProblemE ProblemF ProblemB.txt ProblemC.txt ProblemD.txt ProblemE.txt ProblemF.txt
