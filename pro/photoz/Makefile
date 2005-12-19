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
	- cp *.pro $(IDLSPEC2D_DIR)/pro/photoz
	- cp Makefile $(IDLSPEC2D_DIR)/pro/photoz
	- cp TODO $(IDLSPEC2D_DIR)/pro/photoz

clean :
	- /bin/rm -f *~ core

