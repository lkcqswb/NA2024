# Makefile

# 编译器和编译选项
CXX = g++
CXXFLAGS =

# 目录
SRC_DIR = src
PROBLEMS_DIR = $(SRC_DIR)/problems_in_chapter3
BONUS_DIR = $(SRC_DIR)/bonus

# 可执行文件
PROBLEMS_EXEC = src/problems_1_7/problem1/main src/problems_1_7/problem2/main src/problems_1_7/problem3/main src/problems_1_7/problem4/main src/problems_1_7/problem5/main src/problems_1_7/problem6/main src/problems_1_7/problem7/main
CHAPTER3_EXEC = ProblemA ProblemC ProblemD ProblemE ProblemF
BONUS_EXEC = convergencerate

# 默认目标
all: pardo report

# 编译并运行每个问题目录
pardo: $(PROBLEMS_EXEC) $(CHAPTER3_EXEC) $(BONUS_EXEC)

# 编译每个问题的可执行文件
src/problems_1_7/problem1/main: $(SRC_DIR)/problems_1_7/problem1/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< 
	cd src/problems_1_7/problem1 && ./main

src/problems_1_7/problem2/main: $(SRC_DIR)/problems_1_7/problem2/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< 
	cd src/problems_1_7/problem2 && ./main

src/problems_1_7/problem3/main: $(SRC_DIR)/problems_1_7/problem3/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< 
	cd src/problems_1_7/problem3 && ./main

src/problems_1_7/problem4/main: $(SRC_DIR)/problems_1_7/problem4/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< 
	cd src/problems_1_7/problem4 && ./main

src/problems_1_7/problem5/main: $(SRC_DIR)/problems_1_7/problem5/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< 
	cd src/problems_1_7/problem5 && ./main

src/problems_1_7/problem6/main: $(SRC_DIR)/problems_1_7/problem6/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
	cd src/problems_1_7/problem6 && ./main

src/problems_1_7/problem7/main: $(SRC_DIR)/problems_1_7/problem7/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
	cd src/problems_1_7/problem7 && ./main

# 编译问题章节的可执行文件
ProblemA: $(PROBLEMS_DIR)/ProblemA.cpp
	$(CXX) $(CXXFLAGS) -o src/problems_in_chapter3/ProblemA $< && cd src/problems_in_chapter3 && ./ProblemA

ProblemC: $(PROBLEMS_DIR)/ProblemC.cpp
	$(CXX) $(CXXFLAGS) -o src/problems_in_chapter3/ProblemC $< && cd src/problems_in_chapter3 && ./ProblemC

ProblemD: $(PROBLEMS_DIR)/ProblemD.cpp
	$(CXX) $(CXXFLAGS) -o src/problems_in_chapter3/ProblemD $< && cd src/problems_in_chapter3 && ./ProblemD

ProblemE: $(PROBLEMS_DIR)/ProblemE.cpp
	$(CXX) $(CXXFLAGS) -o src/problems_in_chapter3/ProblemE $< && cd src/problems_in_chapter3 && ./ProblemE

ProblemF: $(PROBLEMS_DIR)/ProblemF.cpp
	$(CXX) $(CXXFLAGS) -o src/problems_in_chapter3/ProblemF $< && cd src/problems_in_chapter3 && ./ProblemF

# 编译 bonus
convergencerate: $(BONUS_DIR)/convergencerate.cpp
	$(CXX) $(CXXFLAGS) -o src/bonus/convergencerate $< && cd src/bonus &&./convergencerate

# 编译文档
report:
	cd doc && pdflatex report.tex
	cd doc && bibtex report
	cd doc && pdflatex report.tex
	cd doc && pdflatex report.tex

# 清理目标
clean:
	rm -f src/problems_1_7/problem1/main src/problems_1_7/problem2/main src/problems_1_7/problem3/main src/problems_1_7/problem4/main src/problems_1_7/problem5/main src/problems_1_7/problem6/main src/problems_1_7/problem7/main
	rm -f src/problems_1_7/problem1/*.txt src/problems_1_7/problem1/*.png src/problems_1_7/problem2/*.txt src/problems_1_7/problem2/*.png src/problems_1_7/problem3/*.txt src/problems_1_7/problem3/*.png src/problems_1_7/problem4/*.txt src/problems_1_7/problem4/*.png src/problems_1_7/problem5/*.txt src/problems_1_7/problem5/*.png src/problems_1_7/problem6/*.txt src/problems_1_7/problem6/*.png src/problems_1_7/problem7/*.txt src/problems_1_7/problem7/*.png
	rm -f src/problems_1_7/problem1/main.exe src/problems_1_7/problem2/main.exe src/problems_1_7/problem3/main.exe src/problems_1_7/problem4/main.exe src/problems_1_7/problem5/main.exe src/problems_1_7/problem6/main.exe src/problems_1_7/problem7/main.exe
	rm -f $(PROBLEMS_DIR)/*.txt $(PROBLEMS_DIR)/*.png $(PROBLEMS_DIR)/*.exe $(PROBLEMS_DIR)/*.csv $(BONUS_DIR)/convergencerate
	rm -f $(PROBLEMS_DIR)/ProblemA $(PROBLEMS_DIR)/ProblemC $(PROBLEMS_DIR)/ProblemD $(PROBLEMS_DIR)/ProblemE $(PROBLEMS_DIR)/ProblemF
	rm -f doc/*.bbl doc/*.blg doc/*.aux doc/*.log doc/*.fls doc/*.gz doc/*.fdb_latexmk src/bonus/convergencerate src/bonus/convergencerate.exe

.PHONY: all pardo report clean