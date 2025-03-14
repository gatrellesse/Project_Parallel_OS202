include Make_linux.inc
include Make_msys2.inc
include Make_osx.inc

# Use mpicxx (MPI compiler) if available, otherwise fallback to g++
CXX = g++

CXXFLAGS = -std=c++17 -fopenmp

LIB = -lmingw32
SDL2FLAGS = -lSDL2main -lSDL2
MPIFLAGS = -lmsmpi

INCLUDES = -I ./include -L ./lib -I "C:\Program Files\Microsoft MPI\SDK\Include" -L "C:\Program Files\Microsoft MPI\SDK\Lib\x64"

ifdef DEBUG
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG
CXXFLAGS2 = $(CXXFLAGS)
else
CXXFLAGS2 = $(CXXFLAGS) $(INCLUDES)
CXXFLAGS += -O3 -march=native -Wall
endif

ALL= simulation.exe 

default:	help

all: $(ALL)

clean:
	rm *.o && rm *.exe

.cpp.o:
	$(CXX) $(CXXFLAGS2) -c $^ -o $@	$(SDL2FLAGS)

simulation.exe : display.o display.hpp model.o model.hpp simulationMPI.o
	$(CXX) $(CXXFLAGS2) $^ -o $@ $(LIB) $(MPIFLAGS) $(SDL2FLAGS)

help:
	@echo "Available targets : "
	@echo "    all            : compile all executables"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"
