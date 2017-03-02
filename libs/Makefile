# Makefile for building the third-party libraries.

# Define a destination directory (local) if it is not already defined (as it should be)
LOC_DIR ?=$(CURDIR)/../local

.PHONY: prerequisites
.DEFAULT: prerequisites
prerequisites: folders lemon sdsl #snap

.PHONY: clean-prerequisites
clean-prerequisites: clean-lemon clean-sdsl #clean-snap

.PHONY: folders
folders:
	@echo "* Checking folders..." ; \
	[ -d $(LOC_DIR) ] || mkdir -p $(LOC_DIR) ; \
	[ -d $(LOC_DIR)/include ] || mkdir -p $(LOC_DIR)/include ; \
	[ -d $(LOC_DIR)/lib ] || mkdir -p $(LOC_DIR)/lib ; \
	[ -d $(LOC_DIR)/lemon ] || mkdir -p $(LOC_DIR)/lemon ; \
	[ -f $(LOC_DIR)/fasta.o ] || echo "* Fasta"

.PHONY: lemon
lemon: folders $(LOC_DIR)/lemon/include/

LEMON_NAME:=lemon-1.3.1
$(LOC_DIR)/lemon/include/:
	@echo "* Lemon Library" ;	\
	tar -xvf $(LEMON_NAME).tar.gz;	\
	if [ -d $(LEMON_NAME)/ ]; then	\
		cd $(LEMON_NAME)/ &&	\
		mkdir build && cd build; \
		sed '3d' ../CMakeLists.txt >> tmp; \
		mv tmp ../CMakeLists.txt; \
		cmake -DCMAKE_INSTALL_PREFIX=$(LOC_DIR)/lemon/ ..; \
		make;	\
		make install; \
	else \
		echo "Lemon library not found in the current directory" ; \
		exit 1 ; \
	fi

.PHONY: clean-lemon
clean-lemon:
	@echo "* Cleaning lemon library..." ; \
	rm -rf $(LEMON_NAME)/; \
	rm -rf $(LOC_DIR)/lemon/

.PHONY: snap
snap: folders $(LOC_DIR)/include/snap-core/

SNAP_NAME:=Snap-2.3
$(LOC_DIR)/include/snap-core/:
	@echo "* SNAP Library" ;	\
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
	@echo "* SDSL-lite Library" ;	\
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
