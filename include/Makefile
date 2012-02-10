###############################################################################
# Sloan Digital Sky Survey III (SDSS-III) -- 2D spectroscopic reduction code
# Code for product: idlspec2d
#
# S. Burles & D. Schlegel
#
# This Makefile & all Makefiles in this product are GNU make compliant.
# Please help keep them that way.  See
# http://www.gnu.org/software/make/manual/make.html
#
# $Id$
#
###############################################################################
#
# This line helps prevent make from getting confused in the case where you
# have a file named 'clean'.
#
.PHONY : clean
#
# This should compile all code prior to it being installed.
#
all : export.h
#
# Some versions of IDL have non-backward compatible export.h files. So use
# whatever is current.
#
export.h :
	@if test -z "$$IDL_DIR"; then \
		echo "IDL_DIR environment variable is not set -- "; \
		echo "it must point to the top of the IDL product directory, "; \
		echo "   e.g. /usr/local/itt/idl"; \
		exit 1; \
	else \
		echo "Linking in IDL export.h, using IDL_DIR=$$IDL_DIR"; \
		test -r "$$IDL_DIR/external/export.h" || (echo "no valid $$IDL_DIR/external/export.h"; exit 1); \
		ln -s $$IDL_DIR/external/export.h; \
	fi

#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) export.h

