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

SUBDIRS = apo2d fplan inspect plate science spec1d spec2d specdb

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :
	@ rm -rf $(IDLSPEC2D_DIR)/pro
	@ mkdir $(IDLSPEC2D_DIR)/pro
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLSPEC2D_DIR)/pro/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done
	- cp Makefile $(IDLSPEC2D_DIR)/pro
	- cp makefile.old $(IDLSPEC2D_DIR)/pro


clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
