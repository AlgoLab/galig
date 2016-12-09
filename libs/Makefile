# Makefile for building the third-party libraries.

# Define a destination directory (local) if it is not already defined (as it should be)
LOC_DIR ?=$(CURDIR)/../local

.PHONY: prerequisites
.DEFAULT: prerequisites
prerequisites: folders snap

.PHONY: clean-prerequisites
clean-prerequisites: clean-snap

.PHONY: folders
folders:
	@echo "* Checking folders..." ; \
	[ -d $(LOC_DIR) ] || mkdir -p $(LOC_DIR)

.PHONY: snap
snap: folders $(LOC_DIR)/include/snap-core/

SNAP_NAME:=Snap-2.3
$(LOC_DIR)/include/snap-core/:
	@echo "* SNAP Lirary" ;	\
	unzip $(SNAP_NAME).zip;	\
	if [ -d $(SNAP_NAME)/ ]; then	\
		cd $(SNAP_NAME)/ &&	\
		$(MAKE) all;	\
		ln -sf ../../thirdparty/$(SNAP_NAME)/snap-core $(LOC_DIR)/include/;	\
		ln -sf ../../thirdparty/$(SNAP_NAME)/glib-core $(LOC_DIR)/include/;	\
	else \
		echo "SNAP library not found in the current directory" ; \
		exit 1 ; \
	fi

.PHONY: clean-snap
clean-snap:
	@echo "* Cleaning SNAP library..." ; \
	rm -rf $(SNAP_NAME)/
