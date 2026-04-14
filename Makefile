# ===== 根目录 Makefile =====

# 可传入参数：相对路径 CASE_DIR，例如 CASE_DIR=cases/LL_1000
CASE_DIR ?= .
ROOT_DIR = $(CASE_DIR)/../..
SRC := $(ROOT_DIR)/2D-DisDyna.cpp
SRC_OMP := $(ROOT_DIR)/2D-DisDyna_parallel.cpp
INCLUDE_DIR := $(ROOT_DIR)/include
TARGET := $(CASE_DIR)/simulate
TARGET_OMP := $(CASE_DIR)/simulate_omp
CODE_SNAPSHOT := $(CASE_DIR)/code/main.cpp

CXX := g++
CXXFLAGS := -O3 -std=c++17 -I$(INCLUDE_DIR) 

.PHONY: build clean

# ===== 默认目标 =====
build: $(TARGET)

# ===== 编译目标 =====
$(TARGET): $(SRC)
	@echo "[INFO] Compiling $< to $(TARGET)..."
	@mkdir -p $(CASE_DIR)/code
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)
	cp $(SRC) $(CODE_SNAPSHOT)

# ===== 编译目标 =====
$(TARGET_OMP): $(SRC_OMP)
	@echo "[INFO] Compiling $< to $(TARGET_OMP)..."
	@mkdir -p $(CASE_DIR)/code
	$(CXX) $(CXXFLAGS) $(SRC_OMP) -o $(TARGET_OMP) -fopenmp
	cp $(SRC_OMP) $(CODE_SNAPSHOT)

# ===== 清理目标（可选）=====
clean:
	@echo "[INFO] Cleaning $(CASE_DIR)..."
	rm -f $(TARGET)
	rm -rf $(CASE_DIR)/code

empty:
	@echo "clean simluation results"
	rm -f outputKD_*
	rm -f disConfig/evl_*
