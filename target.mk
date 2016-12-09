# Based on the idea/skeleton presented at
# http://mad-scientist.net/make/multi-arch.html

.SUFFIXES:

# Directories.
export BASE_DIR:=$(CURDIR)
export SRC_DIR:=$(BASE_DIR)/src
export LOC_DIR:=$(BASE_DIR)/local
export 3RD_DIR:=$(BASE_DIR)/libs
export BIN_DIR:=$(BASE_DIR)/bin
export OBJ_DIR:=$(BASE_DIR)/$(OBJ_DIR)

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile \
	$(MAKECMDGOALS)

.PHONY: $(OBJ_DIR)
$(OBJ_DIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

Makefile : ;
%.mk :: ;

% :: $(OBJ_DIR) ;



# The clean and depclean targets should stay here since they
# delete the directory which the main Makefile is executed from.

.PHONY: clean
clean:
	@echo '* Cleaning objects and binaries...' ; \
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: depclean
depclean: clean
	@echo '* Cleaning local dir...' ; \
	rm -rf $(LOC_DIR) ; \
	echo '* Cleaning pre-requisites...' ; \
	$(MAKE) -C $(3RD_DIR) clean-prerequisites
