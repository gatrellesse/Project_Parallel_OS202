# include Make_linux.inc
include Make_msys2.inc
#include Make_osx.inc

ifdef __linux__
LIB = -lSDL2main -lSDL2 -lmpi
else
LIB = -lmingw32 -lSDL2main -lSDL2 -lmpi
endif

CXX = mpicxx

CXXFLAGS = -std=c++17 -fopenmp -fopenmp-simd -DOMPI_SKIP_MPICXX
# MPIFLAGS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -L/usr/lib/x86_64-linux-gnu 

# SDL2 flags
# SDL2FLAGS = -lSDL2

ifdef DEBUG
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG
CXXFLAGS2 = $(CXXFLAGS)
else
CXXFLAGS2 = ${CXXFLAGS} -O2 -march=native -Wall -I ./include -L ./lib
CXXFLAGS += -O3 -march=native -Wall
endif

ALL= simulation.exe 

default:	help

all: $(ALL)

clean:
	@rm -fr *.o *.exe *~

.cpp.o:
	$(CXX) $(CXXFLAGS2) -c $^ -o $@	$(LIB)

simulation.exe : display.o display.hpp model.o model.hpp simulation.o
	$(CXX) $(CXXFLAGS2) $^ -o $@ $(LIB)

help:
	@echo "Available targets : "
	@echo "    all            : compile all executables"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS2)"

%.html: %.md
	pandoc -s --toc $< --css=./github-pandoc.css --metadata pagetitle="OS202 - TD1" -o $@
	@echo "    CXXFLAGS :    $(CXXFLAGS)"
