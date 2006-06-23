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
	- cp *.$(SO_EXT) $(IDLSPEC2D_DIR)/lib
	- cp Makefile $(IDLSPEC2D_DIR)/lib
clean :
	- /bin/rm -f *~ core *.so *.a *.dylib

