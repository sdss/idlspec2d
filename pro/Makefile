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

SUBDIRS = spec2d

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :
	@echo "You should be sure to have updated before doing this."
	@echo "You will be installing in \$$IDLSPEC2D_DIR=$$IDLSPEC2D_DIR"
	@echo ""
	@echo "I'll give you 5 seconds to think about it"
	@echo sleep 5
	@echo
	@if [ "$(IDLSPEC2D_DIR)" = "" ]; then \
		echo You have not specified a destination directory >&2; \
		exit 1; \
	fi 
	@ rm -rf $(IDLSPEC2D_DIR)
	@ mkdir $(IDLSPEC2D_DIR)
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLSPEC2D_DIR)/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
