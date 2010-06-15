#!/usr/bin/env python

import os, sys

""" Miscellaneous classes for sos programs """

####
class Consts:
	"""Holds various constants use by the sos programs"""
	
	def __init__(self):
		self.lockFileBase = "sos_runner"
		self.dieFileName  = "sos_die.die.die"
		self.logName      = "sos_log"				
		self.configName   = "sos_config.conf"
		self.MJDGlob      = "[0-9]"*5
		self.versionFile  = "sos_version"
	
	def __repr__(self):
		return self.__str__()
			
	def __str__(self):
		return ("Consts:\n"    +
		        "lockFile:   " + self.lockFileBase + "\n" +
		        "dieFile:    " + self.dieFileName + "\n" +
		        "logName:    " + self.logName + "\n" +
		        "configName: " + self.configName + "\n" +
		        "versionFile: " + self.versionFile + "\n" +
		        "MJDGlob:    " + self.MJDGlob);
		
	

####
class Config:
	"""Holds the configuration information read from the command line or ini"""
	
	def __init__(self):
		self.MJD        = "0"
		self.exposure   = None
		self.fitsDir    = "."
		self.plugDir    = "."
		self.controlDir = "."
		self.logDir     = "."
		self.logLevel   = 40
		self.iname      = ""
		self.globs      = ["*"]
		self.command    = ""
		self.pollDelay  = 60
		self.nosvn      = False
		self.bookkeep   = False
		self.nice       = False
		self.platedb    = False
		self.redo       = False
		
	def __repr__(self):
		return self.__str__()

	def __str__(self):
		return ("Config:\n"    +
		        "MJD:        " + self.MJD + "\n" +
		        "Exposure    " + str(self.exposure) + "\n" +
		        "fitsDir:    " + self.fitsDir + "\n" +
		        "plugDir:    " + self.plugDir + "\n" +
		        "controlDir: " + self.controlDir + "\n" +
		        "logDir:     " + self.logDir + "\n" +
		        "logLevel:   " + str(self.logLevel) + "\n" +
		        "iname:      " + self.iname + "\n" +
		        "glob:       " + str(self.globs) + "\n" +
		        "command:    " + self.command + "\n" +
		        "pollDelay:  " + str(self.pollDelay) + "\n" +
		        "nice:       " + str(self.nice) + "\n" +
		        "plateDb:    " + str(self.platedb) + "\n" +
		        "redo:       " + str(self.redo) + "\n" +
		        "bookkeep:   " + str(self.bookkeep) + "\n" +
		        "NoSvn:      " + str(self.nosvn));
			       

####	
class PollWorker:
	"""Holds information for each poller"""
	
	def __init__(self):
		self.workerNumber = 0
		self.glob       = "*"
		self.fileCount  = 0


	def __repr__(self):
		return self.__str__()

	def __str__(self):
		return ("PollWorker:\n" +
		        "workerNumber: " + str(self.workerNumber) + "\n" +
		        "glob:         " + self.glob + "\n" +
		        "fileCount:    " + str(self.fileCount));
		
