#$Id$
#
#  Copyright (c) 1988-1998, Research Systems Inc.  All rights reserved.
#  This software includes information which is proprietary to and a
#  trade secret of Research Systems, Inc.  It is not to be disclosed
#  to anyone outside of this organization. Reproduction by any means
#  whatsoever is  prohibited without express written permission.
#
#  makefile for CALL_EXTERNAL examples, UNIX platforms
#
#  This is a recursive makefile. The default target "all" detects which
#  operating system you are on, sets system specific macros, and then
#  calls make again to build the target "libs".  Specifying a target
#  other than "all","clean", or "tidy" will NOT set all required macros
#  and will probably produce incorrect results.
#
SHELL=/bin/sh
# 	The following macro might need to be edited if the location
#	of IDL is not the standard location on your system.
IDL_DIR = /usr/local/rsi/idl


# The following macros establish general defaults. They are overridden
# in the "all" rule as needed to satisfy a given platform's requirements.

CC		= cc
C_FLAGS		= -I$(IDL_DIR)/external -c $(CFLAGS)
LD 		= ld
SHELL		= /bin/sh
#platform specific C compiler flags (redefined below)
X_CFLAGS	=
#platform specific linker flags
X_LD_FLAGS	=

SO_EXT		=so

.c.o :
	$(CC) $(C_FLAGS) $(X_CFLAGS) $*.c


# The following is the default entry point. This section will determine
# what system we are on, set the correct flags and call this same makefile
# again with the correct flags.

all :
	@echo "OS type detected: "`uname`
	@case `uname` in \
	"SunOS") if [  -d /proc ]; then \
			make libspec2d.$(SO_EXT) \
			"SO_EXT=so" \
			"CC=cc" \
			"LD=cc" \
			"X_CFLAGS=-K pic " \
			"X_LD_FLAGS=-G "; \
		fi \
	;; \
	"HP-UX") make libspec2d.$(SO_EXT) \
		"SO_EXT=sl" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=+z -Aa -D_HPUX_SOURCE" \
		"X_LD_FLAGS=-b";; \
	"IRIX" ) make libspec2d.$(SO_EXT) \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=-32 -KPIC" \
		"X_LD_FLAGS=-shared -32";; \
	"IRIX64" ) make libspec2d.$(SO_EXT) \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=-32 -KPIC" \
		"X_LD_FLAGS=-shared -32";; \
	"OSF1" )  make libspec2d.$(SO_EXT) \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_LD_FLAGS=-S -shared";; \
	"Linux" ) make libspec2d.$(SO_EXT) \
		"SO_EXT=so" \
		"CC=gcc" \
		"LD=ld" \
		"X_CFLAGS=-fPIC -g -Wall" \
		"X_LD_FLAGS=-shared" ;;\
	*) echo "This system is not supported" ;; \
	esac

objs = extract_boxcar.o extract_profile.o extract_row.o \
        sshiftvec.o trace_crude.o trace_fweight.o evilmath.o

libspec2d.$(SO_EXT): $(objs)
	$(LD) $(X_LD_FLAGS) -o ../lib/libspec2d.$(SO_EXT) $(objs)
	nm -s ../lib/libspec2d.$(SO_EXT)

tidy :
	rm -f *.o

clean : tidy
	rm -f *.so *.sl *.a
