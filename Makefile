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
# Change TEMPLATE_DIR here!
#
INSTALL_DIR = $(IDLSPEC2D_DIR)
#
# Use this shell to interpret shell commands, & pass its value to sub-make
#
export SHELL = /bin/sh
#
# This is like doing 'make -w' on the command line.  This tells make to
# print the directory it is in.
#
MAKEFLAGS = w
#
# This is a list of subdirectories that make should descend into.  Makefiles
# in these subdirectories should also understand 'make all' & 'make clean'.
# This list can be empty, but should still be defined.
#
SUBDIRS = include src
#
# This line helps prevent make from getting confused in the case where you
# have a file named 'clean'.
#
.PHONY : clean doc
#
# This should compile all code prior to it being installed
#
all :
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f all ; done
#
# Used to (re)make the documentation files
#
doc :
	$(MAKE) -C doc clean
	$(MAKE) -C doc all
#
# Install things in their proper places in $(INSTALL_DIR)
#
install :
	@ echo "You should be sure to have updated before doing this."
	@ echo ""
	@ if [ "$(INSTALL_DIR)" = "" ]; then \
		echo You have not specified a destination directory >&2; \
		exit 1; \
	fi
	@ if [ -e $(INSTALL_DIR) ]; then \
		echo The destination directory already exists >&2; \
		exit 1; \
	fi
	@ echo ""
	@ echo "You will be installing in \$$INSTALL_DIR=$(INSTALL_DIR)"
	@ echo "I'll give you 5 seconds to think about it"
	@ sleep 5
	@ echo ""
	@ rm -rf $(INSTALL_DIR)
	@ mkdir $(INSTALL_DIR)
	@ cp -Rf . $(INSTALL_DIR)

#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) *~ core
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f clean ; done
