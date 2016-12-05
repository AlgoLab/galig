CC=g++
CFLAGS=-std=c++11
SNAP_FLAGS=-fopenmp ./libs/Snap/snap-core/Snap.o -I ./libs/Snap/snap-core/ -I ./libs/Snap/glib-core/ -lrt
BIN=./bin/
SRC=./src/

main: $(SRC)main.cpp
	$(CC) $(CFLAGS) $(SNAP_FLAGS) -o $(BIN)main $(SRC)main.cpp
