###############################################################################
# Sloan Digital Sky Survey (SDSS) -- 2D spectroscopic reduction code
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh


INC = ../include
LIB = $(IDLSPEC2D_DIR)/lib
CFLAGS  = -I$(IDL_DIR)/external $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC)
#
#.c.o :
%.o : %.c
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $<
#
OBJECTS = \
	extract_boxcar.o \
	extract_row.o \
	sshiftvec.o

all : $(LIB)/libspec2d.$(SO_EXT)


$(LIB)/libspec2d.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libspec2d.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libspec2d.$(SO_EXT)

#
# Install things in their proper places in $(IDLSPEC2D_DIR)
#
install :
	- cp *.c $(IDLSPEC2D_DIR)/src
	- cp Makefile $(IDLSPEC2D_DIR)/src
#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) *~ core *.o
