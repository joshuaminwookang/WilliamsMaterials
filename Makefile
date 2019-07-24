TOPDIR := $(shell pwd)
TXTS := $(shell find $(TOPDIR) -name "*.txt")
FIGS := $(shell find $(TOPDIR) -name "*.fig")

clean :
	@rm -f $(TXTS) $(FIGS) 
	@echo "Removed all .txt .fig .mat files."

