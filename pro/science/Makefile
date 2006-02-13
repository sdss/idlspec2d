###############################################################################
# Sloan Digital Sky Survey (SDSS) -- 2D spectroscopic reduction code
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c
#
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I../include

SUBDIRS = ism kinematics lickindex lya qso

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :
	@ rm -rf $(IDLSPEC2D_DIR)/pro/science
	@ mkdir $(IDLSPEC2D_DIR)/pro/science
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLSPEC2D_DIR)/pro/science/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done
	- cp *.pro $(IDLSPEC2D_DIR)/pro/science
	- cp Makefile $(IDLSPEC2D_DIR)/pro/science
	- cp README $(IDLSPEC2D_DIR)/pro/science

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
