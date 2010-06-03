#!/usr/bin/env python

import os, sys
import sxpar, putils

""" 
filecheck.py:

filecheck.py checks some attribute of a file and returns "true" or "false" based
on if the attribute is...  true or false.

This file it pretty much a grab bag of hacks.  To do these tests correctly would
either require adding dependencies just for these simple tests or duplicating lots
of code.  

The input files can be either normal or gz files.

Written by Gary Kushner (LBL).  Dec 2009.

"""

####
##   Globals Block. 
####

#

####
def usage():
	"""Display usage and exit"""
	
	usageCMD = os.path.basename(sys.argv[0])

	print usageCMD + "\n\nusage: science|boss|test file"
	print \
	"""
	science:
	  return "true" if the fits file is a science frame.  This
	  is determined by flavor=science in the header.  If flavor
	  is not in the header, "false" is returned.
	
	test:
	  return "true" if the fits file is a test frame.  This is
	  determined by quality=test in the header.  If quality
	  is not in the header, "false" is returned
	
	excellent:
	  return "true" if the fits file is a excellent frame.  This is
	  determined by quality=excellent in the header.  If quality
	  is not in the header, "true" is returned
	
	boss:
	  return "true" if the plPlugMapM file is a boss frame.  
	  this is determined by instrument=boss in the header.
	  If instrument is not in the header, "false" is returned.
	
	File can be uncompressed or gz.
	"""
	sys.exit(1)


####
def science(fits):
	"""return True if the fits file is a science frame"""
	v = sxpar.sxpar(fits, "flavor")
	if len(v) == 0:
		return False
	return v[0].lower() == "science"
	


####
def excellent(fits):
	"""return True if the fits file is an excellent frame"""
	v = sxpar.sxpar(fits, "quality")
	if len(v) == 0:
		return True
	return v[0].lower() == "excellent"
	
	
####
def boss(fits):
	"""return True if the plugmap file is for a boss frame"""
	
	key   = "instruments"
	keyl  = len(key)
	value = "boss"
	
	lines = putils.openRead(fits).readlines();
	for line in lines:
		if line[0:keyl].lower() == key:
			return line.lower().count(value) > 0
	return False
		
####
def test(fits):
	"""return True if the fits file is a test frame"""

	v = sxpar.sxpar(fits, "quality")
	if len(v) == 0:
		return False
	return v[0].lower() == "test"


####
def main(argv):
	"""Check a fits file and display 'true' or 'false'"""
	
	if len(argv) != 2:
		usage()
	cmd  = argv[0]
	fits = argv[1]

	if cmd == "science":
		rv = science(fits)
	elif cmd == "test":
		rv = test(fits)
	elif cmd == "excellent":
		rv = excellent(fits)
	elif cmd == "boss":
		rv = boss(fits)
	else:
		usage()

	if rv == True:
		print "true"
	else:
		print "false"
		
### Start of script

if __name__=='__main__':
	main(sys.argv[1:])

