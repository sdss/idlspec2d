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
	- cp *.pro $(IDLSPEC2D_DIR)/pro/science/kinematics
	- cp sp* $(IDLSPEC2D_DIR)/pro/science/kinematics
	- cp Makefile $(IDLSPEC2D_DIR)/pro/science/kinematics

clean :
	- /bin/rm -f *~ core

