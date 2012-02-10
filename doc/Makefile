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
all : www/idlspec2d_doc.html

#
# Replace mk_doc
#
www/idlspec2d_doc.html :
	idl -e 'make_html_help, getenv("IDLSPEC2D_DIR")+"/pro/*/", "www/idlspec2d_doc.html", /strict, title="IDL Help for IDLSPEC2D", /verbose, version="idlspec2d"'

#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) www/idlspec2d_doc.html

