# Makefile

# 编译器和编译选项
CXX = g++
CXXFLAGS = -O2 -Wall

# 目录
SRC_DIR = src
PROBLEMS_DIR = $(SRC_DIR)/problems_in_chapter3
BONUS_DIR = $(SRC_DIR)/bonus

# 可执行文件
PROBLEMS_EXEC = problem1/main problem2/main problem3/main problem4/main problem5/main problem6/main problem7/main
CHAPTER3_EXEC = ProblemA ProblemC ProblemD ProblemE ProblemF
BONUS_EXEC = convergencerate

# 默认目标
all: pardo report

# 编译并运行每个问题目录
pardo: $(PROBLEMS_EXEC) $(CHAPTER3_EXEC) $(BONUS_EXEC)

# 编译每个问题的可执行文件
$(PROBLEMS_EXEC): %: $(SRC_DIR)/problems_in_chapter3/problems%/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< && ./$@

$(CHAPTER3_EXEC): %: $(PROBLEMS_DIR)/Problem%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< && ./$(CHAPTER3_EXEC)

# 编译 bonus
$(BONUS_EXEC): $(BONUS_DIR)/convergencerate.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< && ./$(BONUS_EXEC)

# 编译文档
report:
	cd doc && pdflatex report.tex
	cd doc && @bibtex report
	cd doc && @pdflatex report.tex
	cd doc && @pdflatex report.tex

# 清理目标
clean:
	rm -f $(SRC_DIR)/*/main $(SRC_DIR)/problems_in_chapter3/* $(BONUS_DIR)/convergencerate
	rm -f $(SRC_DIR)/*/*.o $(SRC_DIR)/*/*.txt $(SRC_DIR)/*/*.png
	rm -f doc/*.aux doc/*.log doc/*.pdf

.PHONY: all pardo report clean