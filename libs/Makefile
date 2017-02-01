# Makefile for building the third-party libraries.

# Define a destination directory (local) if it is not already defined (as it should be)
3RD_DIR ?=$(CURDIR)/../libs

.PHONY: prerequisites
.DEFAULT: prerequisites
prerequisites: folders backwardMEM snap

.PHONY: clean-prerequisites
clean-prerequisites: clean-pre

.PHONY: folders
folders:
	@echo "* Checking folders..." ; \
	[ -d $(3RD_DIR) ] || mkdir -p $(3RD_DIR) ; \


.PHONY: backwardMEM
backwardMEM: folders $(3RD_DIR)/../backwardMEM/bin/backwardMEM

bMEM_NAME:=backwardMEM
$(3RD_DIR)/../backwardMEM/bin/backwardMEM:
	@echo "* backwardMEM" ;			\
	cd $(3RD_DIR)/../$(bMEM_NAME)/ &&	\
	make prerequisites &&			\
	make					\

.PHONY: snap
snap: folders $(3RD_DIR)/Snap-2.3/snap-core/

SNAP_NAME:=Snap-2.3
$(3RD_DIR)/Snap-2.3/snap-core/:
	@echo "* SNAP Lirary" ;	\
	unzip -o $(SNAP_NAME).zip;	\
	cd $(SNAP_NAME)/ &&	\
	make all;

.PHONY: clean-pre
clean-pre:
	@echo "* Cleaning Prerequisites..." ; \
	cd $(3RD_DIR)/../$(bMEM_NAME)/ &&	\
	make clean && \
	make clean-prerequisites ; \
	rm -rf $(3RD_DIR)/$(SNAP_NAME)
