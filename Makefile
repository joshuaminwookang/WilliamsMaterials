TOPDIR := $(shell pwd)
TXTS := $(shell find $(TOPDIR) -name "*.txt")
FIGS := $(shell find $(TOPDIR) -name "*.fig")
MATS := $(shell find $(TOPDIR) -name "*.mat")	
clean :
	@rm -f $(TXTS) $(FIGS) $(MATS)
	@echo "Removed all .txt .fig .mat files."

