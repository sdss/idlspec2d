###############################################################################
# Sloan Digital Sky Survey (SDSS) -- 2D spectroscopic reduction code
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
all :;
#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :
	- cp * $(IDLSPEC2D_DIR)/bin
	- cp Makefile $(IDLSPEC2D_DIR)/bin

clean :
	- /bin/rm -f *~ core
