# Makefile

# 编译器和编译选项
CXX = g++
CXXFLAGS = -std=c++11 -I./src

# LaTeX 编译器
LATEX = pdflatex

# 源文件和目标文件
SRC_DIR = src
OBJ_DIR = obj
DOC_DIR = doc

# 所有 cpp 文件
CPP_FILES = $(wildcard $(SRC_DIR)/**/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_FILES))

# main.cpp 文件
MAIN_FILES = $(wildcard $(SRC_DIR)/problems_1_7/*/main.cpp)
MAIN_EXECUTABLES = $(patsubst $(SRC_DIR)/problems_1_7/%/main.cpp,$(SRC_DIR)/problems_1_7/%/main,$(MAIN_FILES))

# 头文件
HEADER_FILES = $(wildcard $(SRC_DIR)/ppForm_a_BSpline/*.hpp)

# 文档
REPORT = $(DOC_DIR)/report.tex

# 默认目标
all: $(OBJ_FILES) $(MAIN_EXECUTABLES)

# 编译头文件（只需编译一次）
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 编译每个 main
$(SRC_DIR)/problems_1_7/%/main: $(OBJ_FILES) $(SRC_DIR)/problems_1_7/%/main.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(OBJ_FILES)

# 编译 LaTeX 文档
$(DOC_DIR)/report.pdf: $(REPORT)
	$(LATEX) $<

# 运行可执行文件
run: $(MAIN_EXECUTABLES)
	@for exe in $(MAIN_EXECUTABLES); do \
		echo "Running $$exe:"; \
		$$exe; \
	done

# 清理目标
clean:
	rm -rf $(OBJ_DIR) $(SRC_DIR)/problems_1_7/*/main $(DOC_DIR)/report.pdf

.PHONY: all clean run