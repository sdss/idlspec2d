###############################################################################
# Sloan Digital Sky Survey (SDSS) -- 2D spectroscopic reduction code
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
all :
	idl mk_doc

#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :
	@ rm -rf $(IDLSPEC2D_DIR)/doc/www
	@ mkdir $(IDLSPEC2D_DIR)/doc/www
	- cp www/*.html $(IDLSPEC2D_DIR)/doc/www
	- cp Makefile $(IDLSPEC2D_DIR)/doc

clean :
	- /bin/rm -f *~ core www/*~ www/core

