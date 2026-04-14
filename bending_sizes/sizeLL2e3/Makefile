# Makefile for dislocation_dynamics_refactored.cpp
# include path: ../../include
#
# Usage:
#   make
#   make FS_LIB=1        # if your old GCC/MinGW still needs -lstdc++fs
#   make clean

CXX      ?= g++
TARGET   := dd_sim
SRC      := dislocation_dynamics_refactored.cpp
OBJ      := $(SRC:.cpp=.o)
INCDIR   := ../../include

STD      ?= -std=c++17
OPT      ?= -O2
WARN     ?= -Wall -Wextra -pedantic
OMP      ?= -fopenmp

# Set FS_LIB=1 on old GCC / old MinGW if linker reports undefined references
FS_LIB   ?= 0
ifeq ($(FS_LIB),1)
    FSLINK := -lstdc++fs
else
    FSLINK :=
endif

CXXFLAGS := $(STD) $(OPT) $(WARN) $(OMP) -I$(INCDIR)
LDFLAGS  := $(OMP) $(FSLINK)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

empty:
	rm -rf ./disConfig
	rm -rf ./output
