###############################################################################
# Sloan Digital Sky Survey (SDSS) -- 2D spectroscopic reduction code
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.c
#
CFLAGS  = -I$(IDL_DIR)/external $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC) 

INC = ../include
LIB = $(IDLSPEC2D_DIR)/lib

OBJECTS = \
	evilmath.o \
	extract_boxcar.o \
	extract_profile.o \
	extract_row.o \
	sshiftvec.o \
	trace_crude.o \
	trace_fweight.o

all : $(LIB)/libspec2d.$(SO_EXT)


$(LIB)/libspec2d.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libspec2d.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libspec2d.$(SO_EXT)

#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :;

clean :
	- /bin/rm -f *~ core *.o so_locations
