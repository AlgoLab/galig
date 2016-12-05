CC=g++
CFLAGS=-std=c++11
SNAP_FLAGS=./libs/Snap/snap-core/Snap.o -I./libs/Snap/snap-core/ -I./libs/Snap/glib-core/
SDSL_FLAGS=-I./libs/sdsl/include -L./libs/sdsl/lib
OTHER=-fopenmp -lrt -lsdsl -ldivsufsort -ldivsufsort64
BIN=./bin/
SRC=./src/


main: $(SRC)main.cpp
	$(CC) $(CFLAGS) $(SNAP_FLAGS) $(SDSL_FLAGS) -o $(BIN)main $(SRC)main.cpp $(OTHER)
