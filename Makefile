######
# Makefile (GNU make is assumed!)
#
# Based on the idea/skeleton presented at
# http://mad-scientist.net/make/multi-arch.html
######

export OBJ_DIR := obj
ifneq ($(OBJ_DIR),$(notdir $(CURDIR)))

include target.mk

else
.DEFAULT: all

VPATH = $(SRC_DIR)

INCLUDE_FLAGS:= -I$(LOC_DIR)/include/ \
                -I$(LOC_DIR)/include/snap-core/ \
                -I$(LOC_DIR)/include/glib-core/	\
				-I$(LOC_DIR)/lemon/include/ \
		-I$(BASE_DIR)/include

# Pre-processor flags
CPPFLAGS= $(INCLUDE_FLAGS)

# Common C and C++ flags
CCXXFLAGS:=-std=c++11 -Wall -O3 -DNDEBUG -march=native -Wno-deprecated -ffunction-sections -fdata-sections -fopenmp
#CCXXFLAGS:=-g -std=c++11 -Wall -O0 -march=native -Wno-deprecated -ffunction-sections -fdata-sections -fopenmp
# C-only flags
CFLAGS+= $(CCXXFLAGS)
# C++-only flags
CXXFLAGS+= $(CCXXFLAGS)
# Linker flags
LDFLAGS+=-Wl,--gc-sections -fopenmp

# Define libraries
LIBS:= \
        -L${LOC_DIR}/lib \
		-L${LOC_DIR}/lemon/lib \

######
#
# LIST OF THE PROGRAMS
# For each listed program 'xxx', two variables must be defined:
# - OBJS_xxx = the object files that compose the program
# - LIBS_xxx = the libraries which must be linked
#
PROGRAMS:=main

# analyzegraph SNAP library
OBJS_main = \
	utils.o \
	FastaReader.o \
	bMEM.o \
	SplicingGraph.o \
	MEMsGraph.o \
	main.o

#LIBS_main= $(LIBS) $(LOC_DIR)/include/snap-core/Snap.o -lrt -lsdsl -ldivsufsort -ldivsufsort64 -lemon
LIBS_main= $(LIBS) -lrt -lsdsl -ldivsufsort -ldivsufsort64 -lemon

#
# END List of programs
######

.PHONY: all
all: $(addprefix $(BIN_DIR)/, $(PROGRAMS))

######
#
# Additional dependencies
#$(BIN_DIR)/analyzeseqs: $(BIN_DIR)/_clustalw2

#
# END Additional dependencies
######

# Build the pre-requisites
.PHONY: prerequisites
prerequisites:
	@echo '* Building pre-requisites...' ; \
	$(MAKE) -C $(3RD_DIR) prerequisites

######
#
# Pattern rules.
# 1. build a .o file from a .cpp file with the same name (via CXX)
# 2. build a BIN_DIR/xxx file from OBJS_xxx and LIBS_xxx
#
######

.PRECIOUS: %.o
%.o: %.cpp
	@echo '* CXX $<'; \
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

.SECONDEXPANSION:

$(BIN_DIR)/%: $$(OBJS_%)
	@echo '* LD  $@'; \
        [ -d $(BIN_DIR) ] || mkdir -p "$(BIN_DIR)" ; \
	$(CXX) $(LDFLAGS) -o $@ $(OBJS_$(notdir $@)) $(LIBS_$(notdir $@))


endif
