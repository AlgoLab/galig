CXX=g++
CFLAGS=-std=c++11 -o0 -g
SNAP_FLAGS=./libs/Snap/snap-core/Snap.o -I./libs/Snap/snap-core/ -I./libs/Snap/glib-core/
SDSL_FLAGS=-I./libs/sdsl/include -L./libs/sdsl/lib
OTHER=-fopenmp -lrt -lsdsl -ldivsufsort -ldivsufsort64
BIN=./bin/
SRC=./src/

all: main tuts

main: $(SRC)main.cpp
	$(CXX) $(CFLAGS) $(SNAP_FLAGS) $(SDSL_FLAGS) -o main $(SRC)main.cpp $(OTHER)

tuts: $(SRC)tuts.cpp
	$(CXX) $(CFLAGS) $(SNAP_FLAGS) $(SDSL_FLAGS) -o tuts $(SRC)tuts.cpp $(OTHER)
