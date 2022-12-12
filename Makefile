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

# Common C and C++ flags
CCXXFLAGS:=-Wall -pthread -g
# C-only flags
CFLAGS+= $(CCXXFLAGS)
# C++-only flags
CXXFLAGS+= $(CCXXFLAGS) -std=c++1z
# Linker flags
LDFLAGS+=-pthread

######
#
# LIST OF THE PROGRAMS
# For each listed program 'xxx', two variables must be defined:
# - OBJS_xxx = the object files that compose the program
# - LIBS_xxx = the libraries which must be linked
#
PROGRAMS:=SpliceAwareAligner
OBJS_SpliceAwareAligner = \
	utils.o \
	bMEM.o \
	SplicingGraph.o \
	MEMsGraph.o \
	SpliceAwareAligner.o

LIBS_SpliceAwareAligner= $(LIBS) -lrt -lsdsl -ldivsufsort -ldivsufsort64 -lemon -lz

#
# END List of programs
######

.PHONY: all
all: $(addprefix $(BIN_DIR)/, $(PROGRAMS))

##########################
# 3RD PART PREREQUISITES #
##########################
.PHONY: prerequisites
prerequisites: salmon

.PHONY: clean-prerequisites
clean-prerequisites: clean-salmon

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
	echo '* CXX $<'; \
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

.SECONDEXPANSION:

$(BIN_DIR)/%: $$(OBJS_%)
	echo '* LD  $@'; \
        [ -d $(BIN_DIR) ] || mkdir -p "$(BIN_DIR)" ; \
	$(CXX) $(LDFLAGS) -o $@ $(OBJS_$(notdir $@)) $(LIBS_$(notdir $@))


endif
