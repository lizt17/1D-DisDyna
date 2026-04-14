# Top-level Makefile
# 功能：
#   make / make all : 所有 case 文件夹同时执行 make && ./dd_sim
#   make build      : 所有 case 文件夹同时执行 make
#   make run        : 所有 case 文件夹同时执行 ./dd_sim
#   make clean      : 所有 case 文件夹执行 make clean
#   make empty      : 所有 case 文件夹执行 make empty

SHELL := /bin/bash

# CASE_DIRS := sizeLL2e3 sizeLL8e3 sizeLL32e3 sizeLL64e3
CASE_DIRS := $(sort $(patsubst %/,%,$(dir $(wildcard */Makefile))))
EXE := ./dd_sim

.PHONY: all build run clean empty cases

all:
	@echo "===== Batch build + run in parallel start ====="; \
	pids=""; \
	for d in $(CASE_DIRS); do \
		if [ -d "$$d" ]; then \
			echo "---- launching $$d ----"; \
			( $(MAKE) -C "$$d" && cd "$$d" && $(EXE) ) & \
			pids="$$pids $$!"; \
		fi; \
	done; \
	status=0; \
	for p in $$pids; do \
		wait $$p || status=1; \
	done; \
	echo "===== Batch build + run done ====="; \
	exit $$status

build:
	@echo "===== Batch build in parallel start ====="; \
	pids=""; \
	for d in $(CASE_DIRS); do \
		if [ -d "$$d" ]; then \
			echo "---- building $$d ----"; \
			( $(MAKE) -C "$$d" ) & \
			pids="$$pids $$!"; \
		fi; \
	done; \
	status=0; \
	for p in $$pids; do \
		wait $$p || status=1; \
	done; \
	echo "===== Batch build done ====="; \
	exit $$status

run:
	@echo "===== Batch run in parallel start ====="; \
	pids=""; \
	for d in $(CASE_DIRS); do \
		if [ -d "$$d" ]; then \
			echo "---- running $$d ----"; \
			( cd "$$d" && $(EXE) ) & \
			pids="$$pids $$!"; \
		fi; \
	done; \
	status=0; \
	for p in $$pids; do \
		wait $$p || status=1; \
	done; \
	echo "===== Batch run done ====="; \
	exit $$status

clean:
	@echo "===== Batch clean start ====="; \
	for d in $(CASE_DIRS); do \
		if [ -d "$$d" ]; then \
			echo "---- cleaning $$d ----"; \
			$(MAKE) -C "$$d" clean; \
		fi; \
	done; \
	echo "===== Batch clean done ====="

empty:
	@echo "===== Batch empty start ====="; \
	for d in $(CASE_DIRS); do \
		if [ -d "$$d" ]; then \
			echo "---- emptying $$d ----"; \
			$(MAKE) -C "$$d" empty; \
		fi; \
	done; \
	echo "===== Batch empty done ====="

cases:
	@echo "Detected case folders:"; \
	for d in $(CASE_DIRS); do echo "$$d"; done