###############################################################################
# Sloan Digital Sky Survey (SDSS) -- 2D spectroscopic reduction code
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = apo2d config fluxfix guider inspect photoz plan plate science science/ism science/kinematics science/lickindex science/lya science/qso spec1d spec2d specdb specflat templates testsuite

all :;

#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLSPEC2D_DIR)/pro/$$f; \
		echo Copying $$f; \
		cp $$f/* $(IDLSPEC2D_DIR)/pro/$$f ); \
	done
	-@ cp Makefile $(IDLSPEC2D_DIR)/pro

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; /bin/rm -f core *~ ); \
	done
