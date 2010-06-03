#!/usr/bin/env python

import sys, os, os, subprocess
import string, imp, time, shlex
import gzip

""" 
putils is a set of miscellaneous python tools.

Written by Gary Kushner (LBL).  Nov 2009.  Latest update April 2010.
"""

def searchPath(name, paths):
	"""Search a path for a name (file, direcory, link, etc).  Return the absolute 
       path to the found file or None"""
	for path in paths:
		if os.path.exists(os.path.join(path, name)):
			return os.path.abspath(os.path.join(path, name))
	return None
	
def loadModuleRaw(module):
	"""import a python module using a raw file name (doesn't need to end in .py)"""
	path = searchPath(module, sys.path)
	if path == None:
		raise ImportError("No module named " + module)
	return imp.load_source(module, path)

def runCommand(cmd, echo=False, logCmd=None, prefix="", shell=False):
	"""Run a command with the option to asynchronously display or log output.
	
	   If shell=False, the cmd needs to be a list, but if you pass in a string
	   it will be parsed into a list.
	
	   echo will echo output to stdout.
	
	   logCmd is a function pointer to use to put the output into a log.
	
	   Returns (return code, output)."""
	
	output = ""
	
	#	Handle the command parsing
	if isinstance(cmd, str) and not shell:
		cmd = [c for c in shlex.split(cmd)]

	#	Call the process
	p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, 
	shell=shell)
	
	#	Process output until process dies
	while True:
		l = p.stdout.readline()
		if not l: break
		output += l
		l = l[:-1]	# yea, only safe on unix...
		if echo:
			print prefix + l
		if logCmd != None:
			logCmd(prefix + l)
	
	return (p.wait(), output)
	
def openRead(filename, mode = "r"):
	"""Open a gzip or normal file for text reading.  Valid modes are 'r' and 'rb'"""
	
	gzSig = '\x1f\x8b'
	
	if mode != 'r' and mode != 'rb':
		raise ValueError("Illegal mode: " + mode)
	
	f = open(filename, mode)
	
	try:
		if (f.read(2) == gzSig):
			f = gzip.open(filename, mode)
	finally:
		f.seek(0)
		
	return f
		
		
		
				
