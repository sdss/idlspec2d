#!/usr/bin/env python

import os, sys, fcntl, time, commands
import logging, logging.handlers, getopt, glob, random
import sos_classes, sxpar, putils
#import fb_classes

""" 
sos_runnerd.py:

sos_runnerd polls a directory structure looking for new files matching a glob and running a command on those
matching files. The directory structure to look for files is root/$MJD.  By default, when run sos_runnerd goes to
the latest MDJ and starts looking for new files. 

Many globs can be specified, and any file matching any specified glob will be processed.  Many instances of
sos_runnerd can be running, for example one per camera.  It is important that only instance be looking for a 
particular kind of file or else files will be processed more than once, things will get overwritten, your
hair will fall out and all sorts of bad things will happen.

sos_runnerd is smart enough not to let two instances with EXACTLY the same globs, specified in EXACTLY the
same order run concurrently.  If you try to start a new instance with exactly the same globs in order then the
process will exit(0).  The only real use of this is to allow chron to start the process as a keep-alive.  The 
lock file to implement this is stored in the control directory.

Everytime the sos_runnerd looks for a new file, if it doesn't find one, it will also look for a later MJD.  If it 
finds one, it will start looking there for new files.  sos_runnerd never looks at previous MJDs and will only look
at one MJD at a time.

sos_runnerd look at command line arguments and a configuration file named sos_config.ini.  Command line arguments
always have precedence.  Parameters and sos_config.ini names are defined in usage().  Configuration files are
not currently supported and probably never will be unless someone asks.

Notes to run:
	In order to run, the environment must be setup:
		setup platedb
		ssh-agent a key that can commit svn
		
	setting up platedb is optional.  runnerd will try and do it automatically, but it
	requires a late version of eups to work and the correct version of platedb must be
	setup as current.

Written by Gary Kushner (LBL).  Oct 2009.

"""



####
def usage():
	"""Display usage and exit"""
	
	usageCMD = os.path.basename(sys.argv[0])

	print """
Parameters and sos_config.ini names are:


iname (-i) : defaut "" : instance name.  Appended to logfiles and other things that need
                         to not clash.  Use the camera name.
glob (-g) : default * : glob to look for.  Can have many.  Each glob will preserve order 
                        processing of files
nosvn (-x) : default not set : Run without doing any svn processing
bookkeep (-k) : default not set : Run bookkeeping every 60 minutes (svn up & ci)
command (-c) : default echo-to-log : Command to run on found files.  
pollDelay (-d) : default 60 seconds : Seconds to wait before interrogating directing
logDir (-l) : default . : Place to place log files
logLevel (-v) : default 60 : 40 = ERROR; -v = 30 = WARNING; -v -v = 20 = INFO; -v -v -v = 10 = DEBUG
controlDir (-o) : default . : Place to find config file and place lock file
fitsDir (-r) : default . : Where to look for new files in MJD subdirectories
plugDir (-p) : default . : Where to look for and put plugmap files in MJD subdirs
nice (-n) : default False : Run commands with nice
platedb (-b) : default False : Exit if can't setup platedb
initialMJD (-m) : default latest : MJD to start looking for new files. Note:  With the -e option,
                                   this option will process files and in the specified MJD and then
                                   go to the newest MJD, NOT the next one.
redo (-e) : default not set : Start with existing the files in the latest or given MJD.  Normally only
                              new files are processed.  Usefull for starting runnerd after observing 
                              has started, or along with -m for a poor man's redo.

For the command, the following substitutions can be made:
   %%f   for the globbed (fits) file name w/o path information
   %%qf  for the fully qualified globbed (fits) file name.
   %%pf  for the path to the globbed (fits) file.
   %%p   for the plugmap file w/o path information
   %%qp  for the fully qualified plugmap file name
   %%pp  for the path to the plugmap file.
   %%m   for the current mjd
	"""
	sys.exit(1)


####
def screamAndDie(msg):
	"""Log a message and then exit"""

	log = logging.getLogger(sos_classes.Consts().logName)
	log.critical(msg)
	log.critical("GOODBYE!")
	
	#	These lines may cause an exception if running headless.  But at this point... so what.
	print >> sys.stderr, msg
	print >> sys.stderr, "GOODBYE!"
	sys.exit(1)

####
def globsToString(cfg):
	"""Return a string representation of all the globs
	
	Currently the string is just a concatinated list of the globs give.  This gives lots of
	characters that are fine for Linux and Darwin but not so much every OS.  This could easily
	be changed to a hash code or translated string if needed.
	"""
	
	s = ""
	for g in cfg.globs:
		s += g
		
	return s
	
	
####
def oneInstanceCheck(cfg, log):
	"""Only one instance of this daemon should be running.  The normal setup is to have chron try
	   and start it every so often and if it is still running, then this procedure will abort.  It
	   writes a message to stdout, but that usually should be shipped to >dev/null"""
	
	lockFile = os.path.join(cfg.controlDir, sos_classes.Consts().lockFileBase)
	lockFile += "-" + globsToString(cfg)
	lockFile += ".lock"
	
	log.info("Lock file is " + lockFile)
	
	lock = open(lockFile, 'w')
	try:
		fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
	except IOError as (errno, errstr):
		log.info("oneInstanceCheck failed gracefully.")
		print >> sys.stderr, "oneInstanceCheck failed gracefully."
		sys.exit(0)
		
	return lock
 

####
def writeVersionInfo(cfg, log):
	"""Write a version string to a file"""
	
	verFile = os.path.join(cfg.controlDir, sos_classes.Consts().versionFile)
	rc = commands.getstatusoutput("idlspec2d_version")
	f = open(verFile, "w")
	f.write(time.ctime() + " " + rc[1] + "\n")
	f.close()
	log.info("Version is %s" % rc[1])
	
	

####
def isitTimeToDie(cfg, log):
	"""check dirFileName to see if it is time to gracefully give up the ghost"""
	
	dieFile = os.path.join(cfg.controlDir, sos_classes.Consts().dieFileName)
	
	if os.path.isfile(dieFile):
		log.critical("Time to die.")
		log.critical("GOODBYE!")
		sys.exit(0)

####
def setupPlateDbAndIdlspec2d(cfg):
	"""Setup plateDb.  By Default do nothing on failure because older versions of eups don't work this way."""

	try:
		import eups
		Eups = eups.Eups(verbose=0)
		cmds = eups.setup(Eups, "platedb", eups.Current())
		if 'false' in cmds and cfg.platedb:
			screamAndDie("Could not setup platedb")
# idlspec2d needs to already be setup or we wouldn't be running...
#		cmds = eups.setup(Eups, "idlspec2d", eups.Current())
#		if 'false' in cmds and cfg.platedb:
#			screamAndDie("Could not setup idlspec2d")
	except:
		if cfg.platedb:
			screamAndDie("Could not setup platedb and/or idlspec2d (except)")

####
def parseConfigFile(cfg):
	"""Parse the config file"""
	

####
def parseCmdLine(cfg):
	"""Parse command line arguments"""
	
	globs    = []		# fill in the command line globs
	verbose  = 0
	
	# parse with options
	try:
		opts, pargs = getopt.gnu_getopt(sys.argv[1:], "i:g:c:d:l:vc:p:r:m:z:o:xnbek")
	except Exception as e:
		print "Illegal option specified."
		print " "
		print str(e)
		print " "
		usage()
	
	if len(pargs) != 0:
		print "All arguments should be parameters (start with '-')"
		print "found " + str(pargs)
		usage()
	
	#	Fill in the config
	for (opt, value) in opts:
		if opt == "-i":
			cfg.iname = value
		if opt == "-g":
			globs.append(value)
		if opt == "-c":
			cfg.command = value
		if opt == "-d":
			cfg.pollDelay= int(value)
		if opt == "-l":
			cfg.logDir = value
		if opt == "-v":
			verbose += 1
		if opt == "-o":
			cfg.controlDir = value
		if opt == "-r":
			cfg.fitsDir = value
		if opt == "-p":
			cfg.plugDir = value
		if opt == "-m":
			cfg.MJD = value;
		if opt == "-z":
			cfg.exposure = value
		if opt == "-x":
			cfg.nosvn = True
		if opt == "-k":
			cfg.bookkeep = True
		if opt == "-n":
			cfg.nice = True
		if opt == "-b":
			cfg.platedb = True
		if opt == "-e":
			cfg.redo = True
			
	#	Any globs override default
	if (len(globs) != 0):
		cfg.globs = globs
	#	Don't want to apply -v on each call, so always start with a base
	if (verbose > 0):
		cfg.logLevel = max(1, sos_classes.Config().logLevel - verbose * 10)
		
	#	Display config values on any verbosity
	if (verbose > 1):
		print "Config values: \n" + str(cfg)


####
def initializeParms():
	"""Initialize all the parameters."""
	cfg = sos_classes.Config();
	
	#	Parse command line to get config.ini information
	parseCmdLine(cfg)
	#	Parse config.ini to get new defaults
	parseConfigFile(cfg)
	#	Parse command line again to give command line precedence 
	parseCmdLine(cfg)
	
	return cfg


####
def initializeLogger(cfg):
	"""Startup logging and set the level"""
	
	lname = os.path.join(cfg.logDir, sos_classes.Consts().logName)
	if cfg.iname != "":
		lname += "-" + cfg.iname
	print "Starting to log to " + lname
	
	log = logging.getLogger(sos_classes.Consts().logName)
	h  = logging.handlers.TimedRotatingFileHandler(lname, 'midnight', 1, 5)
	hc = logging.handlers.TimedRotatingFileHandler(lname + "-error", 'midnight', 1, 5)
#	h = logging.handlers.RotatingFileHandler(lname, maxBytes=1024*1024 * 25, backupCount=3)
	f = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s")
	h.setFormatter(f)
	hc.setFormatter(f)
	h.setLevel(cfg.logLevel)
	hc.setLevel(logging.ERROR)
	log.setLevel(cfg.logLevel)
	log.addHandler(h)
	log.addHandler(hc)
	
	log.critical("Hello. " + sys.argv[0] + " started.")
	log.info("Startup Configuration is: \n\n" + str(cfg) + "\n\n")
	
	return log
	
	
####
def createPollWorkers(cfg, log):
	"""Create poll workers"""

	workers = []

	num = 1
	for glob in cfg.globs:
		p = sos_classes.PollWorker()
		p.glob = glob
		p.workerNumber = num
		num += 1
		workers.append(p)
		log.debug("\nnew PollWorker:\n" + str(p))
		
	return workers
		

####
def resetPollWorkers(workers, cfg, log):
	"""Reset the file count on all the workers"""
	
	log.info("Resetting fileCount on all workers.")
	for w in workers:
		w.fileCount = 0
	
####
def initializePollWorkers(workers, cfg, log):
	"""Initialize poll workers with latest file counts"""
	
	for worker in workers:
		worker.fileCount = len(glob.glob(os.path.join(cfg.fitsDir, cfg.MJD, worker.glob)))
		log.debug("\nInitialized PollWorker:\n" +  str(worker))



####
def lsltr(dir, regex="*"):
	"""return a modification-time sorted list of files in dir"""
	
	files = [os.path.join(dir, f) for f in glob.glob(os.path.join(dir,regex))]
	files.sort(key=lambda tm: os.path.getmtime(tm))
	
	return files
	
	
####
def ls(dir, regex="*"):
	"""return a name sorted list of files in dir"""
	
	files = [os.path.join(dir, f) for f in glob.glob(os.path.join(dir,regex))]
	files.sort()
	
	return files
	
	
	

####
def initializeMJD(cfg, log):
	"""Find the correct MJD to start looking for new files.  If the user specifies an MJD just test
	to see if it exists, otherwise, use the latest MJD."""
	
	#	First check for user specified
	if cfg.MJD != "0":
		path = os.path.join(cfg.fitsDir, cfg.MJD)
		if not os.path.isdir(path):
			screamAndDie("Could not find user specified MJD path: " + path)
		log.info("Using user specified MJD " + path)
	else:
		regex = sos_classes.Consts().MJDGlob;
		try:
			log.debug("Looking for initial MJD in " + cfg.fitsDir)
			cfg.MJD = ls(cfg.fitsDir, regex)[-1][-5:]
			log.info("Latest initial MJD found to be " + os.path.join(cfg.fitsDir, cfg.MJD))
		except:
			screamAndDie("Could not find latest MJD in " + cfg.fitsDir)
			

####
def updateMJD(workers, cfg, log):
	"""Check to see if a new MJD exists"""
	
	regex = sos_classes.Consts().MJDGlob;
	try:
		MJD = ls(cfg.fitsDir, regex)[-1][-5:]
		if (MJD == cfg.MJD):
			return
			
		cfg.MJD = MJD[-5:]
		for worker in workers:
			worker.fileCount = 0

		log.info("Latest updated MJD found to be " + os.path.join(cfg.fitsDir, cfg.MJD))
	except:
		screamAndDie("Could not find latest MJD in " + cfg.fitsDir)


####
def svnAdd(uri, cfg, log):
	"""Add a file or dir to svn"""

#	Even if svn processing is turned off, we should add the file!
#	if cfg.nosvn:
#		return
###
		
	log.info("svn adding " + uri)
	rc = commands.getstatusoutput("svn add " + uri)
	log.info(" -> output:\n" + rc[1])
	if rc[0] != 0:
		log.critical("\nCould not add to svn: " + uri + "\n" + rc[1])
	
def svnCommit(uri, cfg, log):
	"""Run commit on a dir"""
	
	if cfg.nosvn:
		return
		
	log.info("svn committing " + uri)
	rc = commands.getstatusoutput("svn commit " + uri + " -m 'committed by sos_runnerd'")
	log.info(" -> output:\n" + rc[1])
	if rc[0] != 0:
		log.critical("\nCommit failed on " + uri + "\n" + rc[1])
		log.critical("Trying to do a cleanup")
		svnCleanup(uri, cfg,log)
		
def svnUp(uri, cfg, log):
	"""Update a dir"""
	
	if cfg.nosvn:
		return
		
	log.info("svn updating " + uri)
	rc = commands.getstatusoutput("svn up " + uri)
	log.info(" -> output:\n" + rc[1])
	if rc[0] != 0:
		log.critical("\nUpdate failed on " + uri + "\n" + rc[1])
	
	
def svnCheck(uri, cfg, log):
	"""Check that we can access the log of the file.  Return False on not able to access."""

	if cfg.nosvn:
		return True
		
	log.info("Checking svn access to " + uri)
	rc = commands.getstatusoutput("svn log " + uri)
	return rc[0] == 0
	 
def svnCleanup(uri, cfg, log):
	"""Do a cleanup on an SVN dir"""

	if cfg.nosvn:
		return
		
	log.info("svn cleanup " + uri)
	rc = commands.getstatusoutput("svn cleanup " + uri)
	log.info(" -> output:\n" + rc[1])
	if rc[0] != 0:
		log.critical("\nCleanup failed on " + uri + "\n" + rc[1])
	
####
def checkPlugMap(file, cfg, log):
	"""
	Get a plugmap file from the database if needed.  Uses the platedb command catPlPlugMapM so
	make sure platedb is setup!
	
	Returns the fully qualified name of the plugmap file
	"""
	
	dirty = False # svn dirty bit
	speclogDir = cfg.plugDir
	plugmapDir = os.path.join(speclogDir, cfg.MJD)

	log.info("Current plugmap directory is " +  plugmapDir)
		
	#	Get plugmap used by file
	try:
		plugmapFullId = sxpar.sxparRetry(file, "NAME", retries = 5)[0]
	except TypeError as t:
		log.critical("\nCould not parse " + file + "\n ->" + str(t))
		return ""

		
	#	Parse plugmap name
	plugmapName   = "plPlugMapM-" + plugmapFullId + ".par"
	plugParse     = plugmapFullId.split("-")
	plugmapId     = plugParse[0]
	plugmapMJD    = plugParse[1]
	plugmapMapId  = str(int(plugParse[2][0:2]))
	plugmapPtg    = "A"
	if len(plugParse[2]) == 3:
		plugmapPtg = plugParse[2][2]
	log.debug(file + " uses plugmap " + plugmapFullId + " with Id " + plugmapId)
	log.debug("  full name of plugmap file is " + plugmapName)
	log.debug("pId=" + plugmapId + ", pMJD=" + plugmapMJD + ", pMapId=" + plugmapMapId + ", pPointing=" + str(plugmapPtg))
	
	#	See if the plugmap $MJD dir exists, if not create it and add it to svn
	if os.path.isdir(plugmapDir):
		log.info("Found existing directory: " + plugmapDir)
	else:
		log.info("Creating " + plugmapDir)
		#	We let an error happen here because there is a race condition where all the runnerd processes 
		#   fight to create the directory.  If the directory really didn't get created, we'll get
		#   a real error below.
		try:
			os.mkdir(plugmapDir)
			# svnAdd(plugmapDir, cfg, log)
			dirty = True
			log.info("Created " + plugmapDir)
		except:
			pass
			    
		
	#	Check if the file exists, if not get it and add it to svn
	plugpath = os.path.join(plugmapDir, plugmapName)
	if os.path.isfile(plugpath):
		log.info("Found existing plugmap file: " + plugpath)
	else:
		log.info("Getting from platdb: " + plugmapName)
		cmd  = "catPlPlugMapM";
		cmd += " -m " + plugmapMJD
		cmd += " -f " + plugmapMapId
		if plugmapPtg != None:
			cmd += " -p " + plugmapPtg
		cmd += " " + plugmapId
		log.debug("Getting plugmap using: " + cmd)
		rc = commands.getstatusoutput(cmd)
		if rc[0] != 0:
			log.critical("Could not get plugmap for Id " + plugmapId + "\nOutput:\n" + rc[1])
		else:
			#	The file will probably get written multiple times, but at least it won't be corrupted 
			flock = open(plugpath, 'a')
			fcntl.flock(flock, fcntl.LOCK_EX)
			f = open(plugpath, "w")
			f.write(rc[1])
			f.close()
			flock.close()
			log.info("Created " + plugpath)
			# svnAdd(plugpath, cfg, log)
			# svnCommit(plugmapDir, cfg, log)  [going to move to only doing commit during bookkeeping]
	
	return os.path.abspath(plugpath)
	

####
def createCMD(fglob, plugPath, cfg):
	"""Create command with substitutions
	
	%%f   for the globbed (fits) file name w/o path information
	%%qf  for the fully qualified globbed (fits) file name.
	%%pf  for the path to the globbed (fits) file.
	%%p   for the plugmap file w/o path information
	%%qp  for the fully qualified plugmap file name
	%%pp  for the path to the plugmap file.
	%%m   for the current MJD
	"""

	qf = os.path.abspath(fglob)
	f  = os.path.basename(qf)
	pf = os.path.dirname(qf)
	
	qp = os.path.abspath(plugPath)
	p  = os.path.basename(qp)
	pp = os.path.dirname(qp)

	cmd = cfg.command

	cmd = cmd.replace("%%f", f)
	cmd = cmd.replace("%%qf", qf)
	cmd = cmd.replace("%%pf", pf)
	cmd = cmd.replace("%%pp", pp)	# this line needs to be before %%p
	cmd = cmd.replace("%%p", p)
	cmd = cmd.replace("%%qp", qp)
	cmd = cmd.replace("%%m", cfg.MJD)
	
	return cmd 


####
def executeCommand(cmd, cfg, log):
	"""Execute command"""
	
	if cfg.nice:
		cmd = "nice " + cmd
		
	log.info("Running: " + cmd)
#	rc = commands.getstatusoutput(cmd)
#	log.info(" -> rc = " + str(rc[0]))
#	log.debug(" -> output:\n" + rc[1])
	prefix = "cmd: "
	try:
		if cfg.nice:
			i = 1
		else:
			i = 0
		prefix = cmd.split()[i] + ": " 
	except:
		pass
		
	(rc, output) = putils.runCommand(cmd, logCmd=log.debug, prefix=prefix)
	log.info(" -> rc = " + str(rc))
#	log.debug(" -> output:\n" + output)
	if rc != 0:
		log.critical("Command Failed(" + str(rc) + "): " + cmd)
		log.critical("Output:\n" + output)


####		
def	processNewBOSSFiles(worker, files, cfg, log):
	"""  Process new fits files
	
	Check to see if the plugmap file exists in the correct location, if it does not then 
	create it (get it from the database).  Then add the appropiate APO command to the
	correctly numbered process list.  
	
	Before the files are processed, they are sorted by name.  We really want the files
	sorted by time, but because of the sequence number, name is the same as time for any
	given camera.
	
	"""
	
	#	Sort files by name to get into the right time order
	files.sort()
	log.info("Sorted file list" + str(files))
	
	for f in files:
		log.info("processing new file: " + f)
		
		#	Pull plugmap from the db if needed
		plugpath = checkPlugMap(f, cfg, log)
		
		#	Create the command and execute it
		cmd = createCMD(f, plugpath, cfg)
#		plname = fb_classes.Consts().processListName
#		plname = os.path.join(cfg.controlDir, plname) + str(worker.workerNumber)
		executeCommand(cmd, cfg, log)


####
def doBookKeeping(cfg, log):
	"""Do SVN Book Keeping.  Command are not retried."""
	
	speclogDir = cfg.plugDir
	
	log.info("Doing a bookkeeping svn update on $SPECLOG_DIR: " + speclogDir)
	svnUp(speclogDir, cfg, log)
	log.info("Doing a bookkeeping svn add on $SPECLOG_DIR: "+ speclogDir)
	svnAdd(os.path.join(speclogDir, cfg.MJD), cfg, log)
	svnAdd(os.path.join(speclogDir, cfg.MJD, "*"), cfg, log)
	log.info("Doing a bookkeeping svn commit on $SPECLOG_DIR: " + speclogDir)
	svnCommit(speclogDir, cfg, log)
		

		
####
def watch(workers, cfg, log):
	"""  Watch for new files
	
	When a new file comes in read the header to look for the plugmap and then check to see
	if the plugmap file exists.  If it doesn't, get the plugmap from the database and put it
	into the proper MJD directory.  Create the proper MJD directory for the plugmap if needed.
	
	Next, check to see if a newer MJD has been created.  If there are no new files and no new 
	MJD then sleep for cfg.pollDelay. 
	
	Note that only the latest MJD is ever checked, so once a new MJD is created only that MJD
	will be checked.  
	"""		
	
	#	We do some book keeping every 60 minutes.  Calculate the number of pauses between book keeping
	#	Add a random factor so all instances don't try and do it at the same time.
	bookKeepingPauses  = 60 * 60 / cfg.pollDelay
#	bookKeepingPauses += int(random.random() * bookKeepingPauses * .2)  # only one bookkeeper for now
	bookKeepingCount   = bookKeepingPauses
	log.info("Setting bookKeepingPauses to " + str(bookKeepingPauses))


	#	Do an initial bookeeping to make svn is up to date before we start
	if cfg.bookkeep:
		log.info("Entering initial bookkeeping mode")
		doBookKeeping(cfg, log)
	  
	while True:
		pause = True
		
		#	Have we been nicely asked to kill ourselves?
		#	Do we have an outstanding die signal 
		isitTimeToDie(cfg, log)

		#	First check for new files
		for worker in workers:
			files = lsltr(os.path.join(cfg.fitsDir, cfg.MJD), worker.glob)
			if len(files) != worker.fileCount:
				pause = False
				new = len(files) - worker.fileCount
				log.info("Found " + str(new) + " new files in " +
				          os.path.join(cfg.fitsDir, cfg.MJD, worker.glob))
				#	File could get deleted...
				if new > 0:
					processNewBOSSFiles(worker, files[-1 * new:], cfg, log)
				worker.fileCount = len(files)
	
		#	Next check for a new MJD.  Don't wait if there's a new MJD
		if updateMJD(workers, cfg, log):  pause = False
		
		#	Pause if asked
		if pause: 
			if cfg.bookkeep:
				bookKeepingCount -= 1
				if bookKeepingCount < 1:
					log.info("Entering book keeping mode")
					bookKeepingCount = bookKeepingPauses
					doBookKeeping(cfg, log)
				else:
					log.info("Bookkeeping in " + str(bookKeepingCount) + " sleeps.")
			log.info("Sleeping for " + str(cfg.pollDelay) + " seconds.")
			time.sleep(cfg.pollDelay)
		

####
def main():
	"""The program"""
	
	#	Make the logger global so we use it for exception handling	
	global logger
	
	config = sos_classes.Config();
	logger = None

	#	A cry for help?
	if len(sys.argv) > 1 and (sys.argv[1] == "-h" or sys.argv[1] == "-?"):
		usage()
		sys.exit(100)

	#	Initialize
	config = initializeParms()
	logger = initializeLogger(config)
	writeVersionInfo(config, logger)
#	setupPlateDbAndIdlspec2d(config)
	
	#	Do we have an outstanding die signal 
	isitTimeToDie(config, logger)
	
	#	Make sure we're the only copy running with our globs.  If we're not, then exit cleanly
	lock = oneInstanceCheck(config, logger)

	#	Check svn access
	if not svnCheck(config.plugDir, config, logger):
		logger.critical("Could not svn access " + config.plugDir)
		print "Could not svn access " + config.plugDir
		print >> sys.stderr, "Could not svn access " + config.plugDir

	#	Find correct MJD to start on
	initializeMJD(config, logger)

	#	Create poll workers and initialize file counts
	pollWorkers = createPollWorkers(config, logger)
	initializePollWorkers(pollWorkers, config, logger)
	if config.redo:
		resetPollWorkers(pollWorkers, config, logger)

	#	Watch for new files.  Forever...  Unless there are exceptions.  Then
	#	try up to 3 times to get it working.  But mostly...  Forever!
	
	crashes = 3
	while crashes > 0:
		try:
			watch(pollWorkers, config, logger)
		except SystemExit:
			raise
		except:
			crashes = crashes - 1
			if crashes > 0:
				logger.exception("!!! Uncaught exception in watch() !!!  Will Retry  !!!")
			else:
				logger.exception("!!! TOO MANY Uncaught exceptions in watch() !!!")
				raise
