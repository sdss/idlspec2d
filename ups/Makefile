###############################################################################
# Sloan Digital Sky Survey (SDSS) -- 2D spectroscopic reduction code
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
default :
	echo "Please invoke this makefile using sdssmake"
all :;

install :
	$(SDSSTOOLS_DIR)/bin/install_ups_table idlspec2d.table \
							$(IDLSPEC2D_DIR)/ups
	- cp Makefile $(IDLSPEC2D_DIR)/ups

clean :
	rm -f *~ core *.bak *.orig *.old .#* #*#
