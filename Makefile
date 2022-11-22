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

INCLUDE_FLAGS:= -I$(BASE_DIR)/lemon/compiled/include/ \
                -I$(BASE_DIR)/sdsl-lite/compiled/include
#-I$(LOC_DIR)/include/ \
#-I$(LOC_DIR)/include/glib-core/	\

# Pre-processor flags
CPPFLAGS= $(INCLUDE_FLAGS)
# Common C and C++ flags
CCXXFLAGS:=-std=c++1z -Wall -O3 -DNDEBUG -march=native -Wno-deprecated -ffunction-sections -fdata-sections -fopenmp
# C-only flags
CFLAGS+= $(CCXXFLAGS)
# C++-only flags
CXXFLAGS+= $(CCXXFLAGS)
# Linker flags
LDFLAGS+=-Wl,--gc-sections -fopenmp
# Define libraries
LIBS:= -L$(BASE_DIR)/lemon/compiled/lib/ \
       -L$(BASE_DIR)/sdsl-lite/compiled/lib

######
#
# LIST OF THE PROGRAMS
# For each listed program 'xxx', two variables must be defined:
# - OBJS_xxx = the object files that compose the program
# - LIBS_xxx = the libraries which must be linked
#
PROGRAMS:=SpliceAwareAlignerParallel
OBJS_SpliceAwareAlignerParallel = \
	utils.o \
	bMEM.o \
	SplicingGraph.o \
	MEMsGraph.o \
	SpliceAwareAlignerParallel.o

LIBS_SpliceAwareAlignerParallel= $(LIBS) -lrt -lsdsl -ldivsufsort -ldivsufsort64 -lemon -lz

#
# END List of programs
######

.PHONY: all
all: $(addprefix $(BIN_DIR)/, $(PROGRAMS))

##########################
# 3RD PART PREREQUISITES #
##########################
.PHONY: prerequisites
prerequisites: lemon sdsl salmon

.PHONY: clean-prerequisites
clean-prerequisites: clean-lemon clean-sdsl clean-salmon

# LEMON ################################################################

.PHONY: lemon
lemon: $(BASE_DIR)/lemon/

LEMON_NAME:=lemon-1.3.1
$(BASE_DIR)/lemon/:
	@echo "* Lemon Library" ; \
	cd $(BASE_DIR) ; \
	wget http://lemon.cs.elte.hu/pub/sources/$(LEMON_NAME).tar.gz ; \
	tar -xvf $(LEMON_NAME).tar.gz ; \
	rm $(LEMON_NAME).tar.gz ; \
	mv $(LEMON_NAME) lemon ; \
	if [ -d lemon/ ] ; then \
	 	cd lemon/ && mkdir build && cd build ; \
		sed '3d' ../CMakeLists.txt > tmp ; \
		mv tmp ../CMakeLists.txt ; \
		cmake -DCMAKE_INSTALL_PREFIX=../compiled/ .. ; \
		make ; \
		make install ; \
	else \
		echo "! lemon folder not found !" ; \
		exit 1 ; \
	fi

.PHONY: clean-lemon
clean-lemon:
	@echo "* Cleaning lemon library..." ; \
	rm -rf $(BASE_DIR)/lemon/

########################################################################

# SDSL-LITE ############################################################

.PHONY: sdsl
sdsl: $(BASE_DIR)/sdsl-lite/compiled

$(BASE_DIR)/sdsl-lite/compiled:
	@echo "* SDSL-lite Library" ; \
	cd $(BASE_DIR)/sdsl-lite ; \
	mkdir compiled ; \
	./install.sh ./compiled

.PHONY: clean-sdsl
clean-sdsl:
	@echo "* Cleaning SDSL library..." ; \
	rm -rf $(BASE_DIR)/sdsl-lite/compiled

########################################################################

# SALMON ###############################################################

.PHONY: salmon
salmon: $(BASE_DIR)/salmon/bin/salmon

SALMON_NAME:=salmon
$(BASE_DIR)/salmon/bin/salmon:
	@echo "* Salmon" ; \
	cd $(BASE_DIR) ; \
	wget https://github.com/COMBINE-lab/salmon/releases/download/v0.12.0/salmon-0.12.0_linux_x86_64.tar.gz ; \
	tar xfz salmon-0.12.0_linux_x86_64.tar.gz ; \
	mv salmon-0.12.0_linux_x86_64 salmon ; \
	rm salmon-0.12.0_linux_x86_64.tar.gz

.PHONY: clean-salmon
clean-salmon:
	@echo "* Cleaning salmon binary..." ; \
	rm -rf $(BASE_DIR)/salmon

########################################################################

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
