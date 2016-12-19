# Makefile for building the third-party libraries.

# Define a destination directory (local) if it is not already defined (as it should be)
LOC_DIR ?=$(CURDIR)/../local

.PHONY: prerequisites
.DEFAULT: prerequisites
prerequisites: folders snap sdsl

.PHONY: clean-prerequisites
clean-prerequisites: clean-snap clean-sdsl

.PHONY: folders
folders:
	@echo "* Checking folders..." ; \
	[ -d $(LOC_DIR) ] || mkdir -p $(LOC_DIR) ; \
	[ -d $(LOC_DIR)/include ] || mkdir -p $(LOC_DIR)/include ; \
	[ -d $(LOC_DIR)/lib ] || mkdir -p $(LOC_DIR)/lib

.PHONY: snap
snap: folders $(LOC_DIR)/include/snap-core/

SNAP_NAME:=Snap-2.3
$(LOC_DIR)/include/snap-core/:
	@echo "* SNAP Lirary" ;	\
	unzip -o $(SNAP_NAME).zip;	\
	if [ -d $(SNAP_NAME)/ ]; then	\
		cd $(SNAP_NAME)/ &&	\
		$(MAKE) all;	\
		ln -sf ../../libs/$(SNAP_NAME)/snap-core $(LOC_DIR)/include/;	\
		ln -sf ../../libs/$(SNAP_NAME)/glib-core $(LOC_DIR)/include/;	\
	else \
		echo "SNAP library not found in the current directory" ; \
		exit 1 ; \
	fi

.PHONY: clean-snap
clean-snap:
	@echo "* Cleaning SNAP library..." ; \
	rm -rf $(SNAP_NAME)/


.PHONY: sdsl
sdsl: folders $(LOC_DIR)/include/divsufsort64.h

SDSL_NAME:=sdsl-lite
$(LOC_DIR)/include/divsufsort64.h:
	@echo "* SDSL-lite Lirary" ;	\
	git clone https://github.com/simongog/sdsl-lite.git;	\
	if [ -d $(SDSL_NAME)/ ]; then	\
		cd $(SDSL_NAME)/ &&	\
		./install.sh $(LOC_DIR)/;	\
	else \
		echo "SDSL library not found in the current directory" ; \
		exit 1 ; \
	fi

.PHONY: clean-sdsl
clean-sdsl:
	@echo "* Cleaning SDSL library..." ; \
	rm -rf $(SDSL_NAME)/
