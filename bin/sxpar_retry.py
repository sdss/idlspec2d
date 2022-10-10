#!/usr/bin/env python

import sys, os, getopt, time
import putils
from sxpar import *
""" 
sxpar:

Simply parse a fits header.  Copied from perl "sxpar by D. Finkbeiner 2001 Dec 20".

Can read uncompressed or gz files.

Written by Gary Kushner (LBL).  Oct 2009.

"""

####
def usage():
	"""Display usage and exit"""
	
	usageCMD = os.path.basename(sys.argv[0])

	print("usage:")
	print("\t%s [-v] fits-file [keyword]" % usageCMD)

	sys.exit(1)

####
def main(argv):
	"""Parse arguments and run the script"""
	
	fitsfile = None
	verbose  = False
	keyword  = None
	
	if len(argv) == 0:
		usage()
		
		
	# parse with options
	try:
		opts, pargs = getopt.gnu_getopt(argv, "v")
	except:
		usage()
	
	if len(pargs) == 0:
		usage()
	
	for (opt, value) in opts:
		if opt == "-v":
			verbose = True
	
	fitsfile = pargs[0]
	if len(pargs) > 1:
		keyword = pargs[1]
		
	output = sxparRetry(fitsfile, keyword, verbose, retry = 60)
	
	for l in output:
		print(l)
	
	
		
if __name__=='__main__':
	main(sys.argv[1:])
		
		
	
