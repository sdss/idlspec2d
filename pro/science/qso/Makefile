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
	- cp *.pro $(IDLSPEC2D_DIR)/pro/science/qso
	- cp sp* $(IDLSPEC2D_DIR)/pro/science/qso
	- cp Makefile $(IDLSPEC2D_DIR)/pro/science/qso

clean :
	- /bin/rm -f *~ core

