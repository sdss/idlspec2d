#!/usr/bin/env python

import os, sys, fcntl

""" Miscellaneous classes for FemtoBatch programs """

####
class Consts:
	"""Holds various constants use by the FemtoBatch programs"""
	
	def __init__(self):
		self.logName         = "fb_log"
		self.configName      = "fb_config.conf"
		self.processListName = "fb_plist"
	
	def __repr__(self):
		return self.__str__()
			
	def __str__(self):
		return ("Consts:\n"    +
		        "logName:    " + self.logName + "\n" +
		        "configName: " + self.configName + "\n" +
		        "plist:      " + self.processListName);
		
	

####
class Config:
	"""Holds the configuration information read from the command line or ini"""
	
	def __init__(self):
		self.controlDir = "."
		self.logDir     = "."
		self.logLevel   = 30
		self.pollDelay  = 60
		self.name       = ""
		self.nice       = False
		
	def __repr__(self):
		return self.__str__()

	def __str__(self):
		return ("Config:\n"    +
		        "controlDir: " + self.controlDir + "\n" +
		        "logDir:     " + self.logDir + "\n" +
		        "logLevel:   " + str(self.logLevel) + "\n" +
		        "pollDelay:  " + str(self.pollDelay) + "\n" +
		        "name:       " + self.name + "\n" +
		        "nice:       " + str(self.nice))
			       

class ProcessList:
	"""Manipulate the process files"""

	def __init__(self, filename):
		self.filename = filename
		
	def append(self, item):
		"""Add an item to the process file"""
		if len(item.rstrip()) < 1: return
		if item[-1] != '\n': item += ('\n')

		f = open(self.filename, 'a')
		fcntl.flock(f, fcntl.LOCK_EX)
		f.write(item)
		f.close

	def pop(self):
		"""Return an item from the process file or None"""
		item = None

		if os.path.isfile(self.filename):
			f = open(self.filename, "r+")
			fcntl.flock(f, fcntl.LOCK_EX)
			items = f.readlines()
			if len(items) > 0:
				item = items[0].rstrip()
				items = items[1:]
				f.seek(0)
				f.truncate()
				f.writelines(items)
			f.close()

		return item
