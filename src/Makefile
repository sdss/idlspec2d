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
# Pass these flags to the C compiler
#
# SDSS_CFLAGS X_CFLAGS X_LD_FLAGS SO_EXT are defined by evilmake
# CCCHK is NOT defined by evilmake but is defined by sdssmake
# CCCHK = -Wall is the standard for Linux
#
INC = $(IDLSPEC2D_DIR)/include
LIB = $(IDLSPEC2D_DIR)/lib
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC)
CCCHK = -Wall
#
# How to make objects from c files
#
%.o : %.c
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $<
#
# Objects to compile
#
OBJECTS = \
	extract_boxcar.o \
	extract_row.o \
	sshiftvec.o

all : $(LIB)/libspec2d.$(SO_EXT)

$(LIB)/libspec2d.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libspec2d.$(SO_EXT) $(OBJECTS)
#	nm -s $(LIB)/libspec2d.$(SO_EXT)

#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) *~ core *.o $(LIB)/libspec2d.$(SO_EXT)
