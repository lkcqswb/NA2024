
CXX = g++


CXXFLAGS = -std=c++11 -Wall


TARGETS = ProblemB ProblemC ProblemD ProblemE ProblemF


report:report.pdf


%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run:$(TARGETS)

report.pdf: report.tex
	@pdflatex $report.tex
	@bibtex main
	@pdflatex $report.tex
	@pdflatex $report.tex


clean:
	rm -f $(TARGETS) report.tex report.aux report.log report.pdf

