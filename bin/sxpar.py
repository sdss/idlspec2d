#!/usr/bin/env python

import sys, os, getopt, time
import putils

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

	print "usage:"
	print "\t%s [-v] fits-file [keyword]" % usageCMD

	sys.exit(1)
	
####
def sxparRetry(fitsfile, keyword = None, verbose = False, retries = 0):
	"""call sxpar with retries.  To see if the fits wasn't finished writing."""
	while True:
		try:
			return sxpar(fitsfile, keyword, verbose)
		except:
			if retries < 0:
				raise
			retries = retries - 1
			time.sleep(1)
			
####
def sxpar(fitsfile, keyword = None, verbose = False):
	"""Parse fits header and return output list"""
	
	isFits = False
	output = []

	if keyword != None:
		keyword = keyword.upper()
	
	f = putils.openRead(fitsfile)
	i = 0
	while True:
		i += 1
		line = f.read(80)
		
		key = line.split("=")[0].strip()
		
		if key == "SIMPLE":
			isFits = True
			continue
		if key == "END":
			break

		if i > 40 and not isFits:
			raise TypeError(fitsfile + " Doesn't look like a fits file -- did not find 'SIMPLE'")

		if keyword == None or key == keyword:
			values  = line.partition("=")[2].partition("/")
#			comment = values[2].strip()
			value   = values[0].strip().strip("'").strip()
			
			if verbose:
				output.append(line)
			elif keyword == None:
				output.append(key + " = " + value)
			else:
				output.append(value)
		#fi
	
	return output


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
		
	output = sxpar(fitsfile, keyword, verbose)
	
	for l in output:
		print l
	
	
		
if __name__=='__main__':
	main(sys.argv[1:])
		
		
	