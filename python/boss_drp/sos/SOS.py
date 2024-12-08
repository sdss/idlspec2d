#!/usr/bin/env python3
import boss_drp
from boss_drp.sos import sos_classes
from boss_drp.utils import putils, sxpar
from boss_drp.utils import sxpar
from boss_drp.utils.hash import create_hash
from boss_drp.prep.readfibermaps import readfibermaps

import argparse
from argparse import ArgumentTypeError
import sys
import logging
import logging.handlers
import subprocess
import time
import os
import re
import glob
import copy
import numpy as np
from astropy.io.fits import getheader
from multiprocessing import Process
import datetime



####
class fs_Config:
    """file sequencer Config Info"""

    def __init__(self, cfg):
        self.fitname  = ""
        self.fitdir   = ""
        self.plugname = ""
        self.plugdir  = ""
        self.run_config = cfg

    def __str__(self):
        return ("fitname:  " + self.fitname + "\n" +
                "fitdir:   " + self.fitdir + "\n" +
                "plugname: " + self.plugname + "\n" +
                "plugdir:  " + self.plugdir + "\n" +
                "------------------------------\n" +
                str(self.run_config));

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
        log.critical("Could not find latest MJD in " + cfg.fitsDir)
        log.critical("GOODBYE!")
        sys.exit(1)

def flavor(cfg):
    """return the flavor of the fits file"""
    fv = sxpar.sxparRetry(os.path.join(cfg.fitdir,cfg.fitname), "flavor", retries = 5)[0].lower()
    if fv == "arc":
        hart = sxpar.sxparRetry(os.path.join(cfg.fitdir,cfg.fitname), "HARTMANN", retries = 5)[0].lower()
        if hart in ['right','left']:
            fv='hart'
    return fv

def plugging(cfg):
    """return the plugging of the fits file"""
    return sxpar.sxparRetry(os.path.join(cfg.fitdir,cfg.fitname), "CONFID", retries = 5)[0].lower()


def previousExposure(cfg, log):
    """return a config for the previous exposure"""

    #    Parse: sdR-b1-00114186.fit.gz
    left  = cfg.fitname[:7]
    right = cfg.fitname[15:]
    exp   = (str(int(cfg.fitname[7:15]) - 1)).zfill(8)

    prevcfg = copy.copy(cfg)
    prevcfg.fitname = left + exp + right
    log.info("previous cfg:\n" + str(prevcfg))
    return prevcfg


def rule1(cfg, log):
    """Handle arc/flat ordering"""
    prevcfg = previousExposure(cfg, log)
    log.info("Exposure Flavor: " + flavor(cfg))
    prvfitsExist = os.path.exists(os.path.join(prevcfg.fitdir, prevcfg.fitname))
    log.info("Previous exposure exists: " + str(prvfitsExist))
    if prvfitsExist:
        log.info("Previous Flavor: " + flavor(prevcfg))
        log.info("Same Plugging: " + str(plugging(cfg) == plugging(prevcfg)))

    #    Handle flats -- process, and arc if was previous
    if flavor(cfg) == "flat":
        log.info("Exposure is a flat")
        processFile(cfg, log,  "flat")
        if prvfitsExist and flavor(prevcfg)  == "arc":
            log.info("Processing previous arc")
            processFile(prevcfg, log,  "arc")
        return True
    return False


####
def processFile(cfg, log, flavor=""):
    """call sos_command on the file.  Will exit with error code if the command failts."""

    cmd  = "sos_command"
    cmd += " -f " + cfg.fitname
    cmd += " -i " + cfg.fitdir
    cmd += " -p " + cfg.plugname
    cmd += " -l " + cfg.plugdir
    cmd += " -s " + cfg.run_config.sosdir
    cmd += " -m " + cfg.run_config.MJD
    if cfg.run_config.fps:
        cmd += " -e "
    if cfg.run_config.nocal:
        cmd += " -a "
    if cfg.run_config.nodb:
        cmd += " -n "
    if cfg.run_config.no_reject:
        cmd += " -r "
    if cfg.run_config.sdssv_sn2:
        cmd += " -v "
    if cfg.run_config.arc2trace:
        cmd += " -t "
    if cfg.run_config.forcea2t:
        cmd += " -o "
    if cfg.run_config.pause:
        cmd += " -j "+sos_classes.Consts().licensePause
    if cfg.run_config.utah:
        cmd += " -c "
        cmd += " -u "
    if cfg.run_config.sn2_15:
        cmd += " -b "
    
    prefix = "sos_command(" + flavor + "): "

    echo = cfg.run_config.termverbose

    i = 0
    while i < 5:
        if i > 0:
            sleep(2)
            log.info("Trying again to get idl license")
        log.info("executing: " + cmd)
        rv = putils.runCommand(cmd, echo=echo, prefix=prefix, logCmd=log.info, limit=10)
        if rv[0] != 0:
            log.info("\nCommand failed with rc = " + str(rv[0]) + "\n")
            sys.exit(1)
        if 'Failed to acquire license.' not in rv[1]: break
        
    test = create_hash(os.path.join(cfg.run_config.sosdir,cfg.run_config.MJD))
    if test:
        log.info("\nsha1sum is locked")




def sos_filesequencer(fitname, fitpath, plugname, plugpath, cfg, log=None):
    """
        Checks if processing a flat if there was an arc before the flat,
            then process the after the flat
        Otherwise process the fits file
    """


    config = fs_Config(cfg)
    config.fitname  = fitname
    config.fitdir   = fitpath
    config.plugname = plugname
    config.plugdir  = plugpath
    #    Rules return true if they processed the file and processing should stop
    #    process rule

    log.info("Checking Rule 1")
    if not rule1(config, log):
        log.info("Passed Rules; let's go!")
        processFile(config, log, flavor(config))

####
def processNewBOSSFiles(worker, files, cfg, log):
    """  Process new fits files

    Get the plugmap name and then add the appropiate sos command to the
    correctly numbered process list.

    Before the files are processed, they are sorted by name.  We really want the files
    sorted by time, but because of the sequence number, name is the same as time for any
    given camera.

    """

    #   Sort files by name to get into the right time order
    files.sort()
    log.info("Sorted file list" + str(files))
    if len(files) == 0:
        log.info('No Files')
    for file in files:
        log.info("\n+++++++++++++++++++++++++++++++++++++++++++ \nprocessing new file: " + file)

        #- Get platetype from header (missing=BOSS)
        try:
            hdr = getheader(file)
        except:
            time.sleep(5)
            hdr = getheader(file)
            

        if 'FLAVOR' not in hdr:
            log.info("Skipping exposure with missing FLAVOR keyword.")
            return
        else:
            flavor = hdr['FLAVOR']

        if 'PLATETYP' in hdr:
            platetype = hdr['PLATETYP'].upper()
        else:
            platetype = 'BHM&MWM'

        #- always process bias and darks, regardless of PLATETYP
        #- for other exposure types, only process BOSS and EBOSS exposures
        if flavor in ('bias', 'dark') or platetype in ('BHM', 'BHM&MWM'):
            plugpath = getPlugMap(file, cfg, log)

            SOS_opts = {'confSummary':plugpath,'ccd':os.path.basename(file).split('-')[1], 'mjd':cfg.MJD,
                        'log':True, 'log_dir': os.path.join(cfg.sosdir,str(cfg.MJD))}
            readfibermaps(topdir=os.path.join(cfg.sosdir,str(cfg.MJD)), SOS=True, SOS_opts=SOS_opts,
                          logger=log, clobber=cfg.clobber_fibermap)

            qf = os.path.abspath(file)
            fitname  = os.path.basename(qf)
            fitpath = os.path.dirname(qf)

            plugPath = os.path.abspath(plugpath)
            plugname  = os.path.basename(plugPath)
            plugpath = os.path.dirname(plugPath)

            sos_filesequencer(fitname, fitpath, plugname, plugpath, cfg, log=log)
        else:
            #- Don't crash if hdr is mangled and doesn't have EXPOSURE
            if 'EXPOSURE' in hdr:
                log.info("Skipping %s exposure %d." % (platetype, hdr['EXPOSURE']))
            else:
                log.info("Skipping %s exposure." % platetype)
            return

def writeVersionInfo(cfg, log):
    """Write a version string to a file"""

    verFile = os.path.join(cfg.controlDir, sos_classes.Consts().versionFile)
    f = open(verFile, "w")
    f.write(time.ctime() + " " + boss_drp.__version__ + "\n")
    f.close()
    log.info("Version is %s" % boss_drp.__version__)
    log.info("IDLSPEC2D Module is %s" % os.getenv('IDLSPEC2D_VER'))

####
def getPlugMap(file, cfg, log):
    """
    Returns the fully qualified name of the plugmap file
    """

    speclogDir = cfg.plugDir


    try:
        obs = os.environ['OBSERVATORY'].lower()
    except:
        log.critical('Not running at APO or LCO, set OBS environmental variable to run')
        return ""
    plugmapDir = os.path.join(speclogDir, obs, 'summary_files')

    #   Get plugmap used by file
    try:
        plugmapFullId = sxpar.sxparRetry(file, "CONFID", retries = 5)[0]
    except TypeError as t:
        log.critical("\nCould not parse " + file + "\n ->" + str(t))
        return ""
    except IndexError:
        log.critical("\nKeyword CONFID not found in " + file)
        return ""
    try:
        if int(plugmapFullId) == -999:
            plugmapFullId = '0'
    except:
        plugmapFullId = '0'
#    try:
#        plugmapDir = os.path.join(plugmapDir, str(int(np.floor(int(plugmapFullId)/100))).zfill(4)+'XX')
#    except:
#        plugmapDir = os.path.join(plugmapDir, str(int(np.floor(int(0)))).zfill(4)+'XX')

    try:
        configgrp = '{:0>3d}XXX'.format(int(plugmapFullId)//1000)
        configdir = '{:0>4d}XX'.format(int(plugmapFullId)//100)
        plugmapDir = os.path.join(plugmapDir, configgrp, configdir)#,'confSummary-'+str(configid)+'.par
    except:
        configgrp = '{:0>3d}XXX'.format(int(0)//1000)
        configdir = '{:0>4d}XX'.format(int(0)//100)
        plugmapDir = os.path.join(plugmapDir, configgrp, configdir)#,'confSummary-'+str(configid)+'.par

#    try:
#        plugmapDir = os.path.join(plugmapDir, str(int(np.floor(int(plugmapFullId)/1000))).zfill(3)+'XXX',
#                                  str(int(np.floor(int(plugmapFullId)/100))).zfill(4)+'XX')
#    except:
#        plugmapDir = os.path.join(plugmapDir, str(int(np.floor(int(0)/1000))).zfill(3)+'XXX',
#                                  str(int(np.floor(int(0)/100))).zfill(4)+'XX')
    log.info("Current confSummary directory is " +  plugmapDir)

    #   Parse plugmap name
    plugmapName   = "confSummaryF-" + plugmapFullId + ".par"
    plugParse     = plugmapFullId.split("-")
    plugmapId     = plugmapFullId #plugParse[1]
    log.debug(file + " uses confSummaryF " + plugmapFullId + " with Id " + plugmapId)
    log.debug("  full name of confSummaryF file is " + plugmapName)
    log.debug("confId=" + plugmapId)# + ", pMJD=" + plugmapMJD)

    #   Check if the file exists, if not get it and add it to svn
    plugpath = os.path.join(plugmapDir, plugmapName)
    if os.path.isfile(plugpath):
        log.info("Found existing confSummaryF file: " + plugpath)
    else:
        plugmapName   = "confSummary-" + plugmapFullId + ".par"
        plugParse     = plugmapFullId.split("-")
        plugmapId     = plugmapFullId #plugParse[1]
        log.info("No confSummaryF file: " + plugpath)
        log.debug(file + " uses confSummary " + plugmapFullId + " with Id " + plugmapId)
        log.debug("  full name of confSummary file is " + plugmapName)
        log.debug("confId=" + plugmapId)# + ", pMJD=" + plugmapMJD)

        plugpath = os.path.join(plugmapDir, plugmapName)

        if os.path.isfile(plugpath):
            log.info("Found existing confSummary file: " + plugpath)
        else: log.critical("Could not get confSummary for Id " + plugmapId)
    return os.path.abspath(plugpath)

def initializePollWorkers(workers, cfg, log):
    """Initialize poll workers with latest file counts"""

    for worker in workers:
        worker.fileCount = len(glob.glob(os.path.join(cfg.fitsDir, cfg.MJD, worker.glob)))
        log.debug("\nInitialized PollWorker:\n" +  str(worker))

def initializeMJD(cfg, log):
    """Find the correct MJD to start looking for new files.  If the user specifies an MJD just test
    to see if it exists, otherwise, use the latest MJD."""

    #   First check for user specified
    if cfg.MJD != "0":
        path = os.path.join(cfg.fitsDir, cfg.MJD)
        if not os.path.isdir(path):
            log.critical("Could not find user specified MJD path: " + path)
            log.critical("GOODBYE!")
        log.info("Using user specified MJD " + path)
    else:
        regex = sos_classes.Consts().MJDGlob;
        try:
            log.debug("Looking for initial MJD in " + cfg.fitsDir)
            cfg.MJD = ls(cfg.fitsDir, regex)[-1][-5:]
            log.info("Latest initial MJD found to be " + os.path.join(cfg.fitsDir, cfg.MJD))
        except:
            log.critical("Could not find the latest MJD path: " + cfg.fitsDir)
            log.critical("GOODBYE!")

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


def initializeLogger(cfg):
    """Startup logging and set the level"""

    lname = os.path.join(cfg.logDir, sos_classes.Consts().logName)
    if cfg.iname != "":
        lname += "-" + cfg.iname

    log = logging.getLogger(sos_classes.Consts().logName)
    rollover = datetime.time(hour=18)
    ext = '.log' if cfg.utah else ''
        
    print("Starting to log to " + lname+ext+" ")
    h  = logging.handlers.TimedRotatingFileHandler(lname+ext,
                                                  when='midnight',
                                                  interval=1, backupCount=5,
                                                  atTime=rollover)
    hc = logging.handlers.TimedRotatingFileHandler(lname + "-error"+ext,
                                                  when='midnight',
                                                  interval=1, backupCount=5,
                                                  atTime=rollover)

    f = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s")
    h.setFormatter(f)
    hc.setFormatter(f)
    h.setLevel(cfg.logLevel)
    hc.setLevel(logging.ERROR)
    log.setLevel(cfg.logLevel)
    log.addHandler(h)
    log.addHandler(hc)

    log.info("Hello. " + sys.argv[0] + " started.")
    log.info("Startup Configuration is: \n\n" + str(cfg) + "\n\n")

    return log

def closeLogger(log):
    handlers = log.handlers[:]
    for handler in handlers:
        log.removeHandler(handler)
        handler.close()

def ls(dir, regex="*"):
    """return a name sorted list of files in dir"""

    files = [os.path.join(dir, f) for f in glob.glob(os.path.join(dir,regex))]
    files.sort()

    return files

def lsltr(dir, regex="*"):
    """return a modification-time sorted list of files in dir"""

    files = [os.path.join(dir, f) for f in glob.glob(os.path.join(dir,regex))]
    files.sort(key=lambda tm: os.path.getmtime(tm))
    files1 = [f for f in files if len(os.path.basename(f).split('.'))==3]

    return files1


def redo(workers, cfg, log):
    """Redo the command for files in the specified MJD"""

    #    Get files
    for worker in workers:
        files = lsltr(os.path.join(cfg.fitsDir, cfg.MJD), worker.glob)
        if cfg.exposure != None:
            allfiles = files
            files = []
            #    should only be one, but I do it this way to be sure things are working.
            #   also, I do the numeric check to avoid worrying about leading zeroes.
            for file in allfiles:
                log.info("Checking exposure number of: " + file)
                exp = re.search("sdR\-..-(\d{8})\.fit.*$", file)
                if exp != None:
                    if int(exp.group(1)) == int(cfg.exposure):
                        log.info("correct exposure number")
                        files.append(file)
        new = len(files)
        log.info("Found " + str(new) + " files in " +
                  os.path.join(cfg.fitsDir, cfg.MJD, worker.glob))
        processNewBOSSFiles(worker, files, cfg, log)


def runner(pollWorkers, config, log):
    """
    Run the check for fits files and processes it
    """
        #   Create poll workers and initialize file counts
    initializePollWorkers(pollWorkers, config, logger)
    #   Watch for new files.  Forever...  Unless there are exceptions.  Then
    #   try up to 5 times to get it working.  But mostly...  Forever!
    crashes = 5
    while crashes > 0:
        try:
            watch(pollWorkers, config, logger)
        except SystemExit:
            raise
        except:
            crashes = crashes - 1
            if crashes > 0:
                logger.exception("!!! Uncaught exception in watch() !!!  Will Retry  !!!")
                time.sleep(5)
            else:
                logger.exception("!!! TOO MANY Uncaught exceptions in watch() !!!")
                raise

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
    tpause = 0
    vpollDelay = 300
    while True:
        pause = True
        #   First check for new files
        for worker in workers:
            files = lsltr(os.path.join(cfg.fitsDir, cfg.MJD), worker.glob)
            if len(files) != worker.fileCount:
                pause = False
                new = len(files) - worker.fileCount
                log.info("Found " + str(new) + " new files in " +
                          os.path.join(cfg.fitsDir, cfg.MJD, worker.glob))
                #   File could get deleted...
                if new > 0:
                    processNewBOSSFiles(worker, files[-1 * new:], cfg, log)
                worker.fileCount = len(files)
                tpause = 0

        #   Next check for a new MJD.  Don't wait if there's a new MJD
        if updateMJD(workers, cfg, log):  pause = False

        #   Pause if asked
        if pause:
            if (tpause % vpollDelay == 0) or (tpause / vpollDelay >= 1):
                log.info(f"Sleeping for {cfg.pollDelay} seconds. (outout every {vpollDelay}s)")
                tpause = 0
            tpause += cfg.pollDelay
            time.sleep(cfg.pollDelay)

def SOS(CCD, exp=None, mjd=None, catchup=False, redoMode=False,systemd=False, nodb=False,
        no_gz=False, no_reject=False, clobber_fibermap=False, sdssv_sn2=False,
        arc2trace=False, forcea2t=False, pause = False, test=False, utah=False,
        termverbose = False, sn2_15=False):
    """
    The SOS controller for both manual runs and systemd tasks
    """
    global logger
    for i, ex in enumerate(exp):
        if i > 1:
            pause = False
        config = sos_classes.Config();
        config.setup(CCD = CCD, mjd=mjd, exp=ex,
                     redo=redoMode, catchup=catchup, test=test, systemd=systemd,
                     no_gz=no_gz, nodb=nodb, no_reject=no_reject, sdssv_sn2=sdssv_sn2,
                     pause=pause, arc2trace=arc2trace, forcea2t=forcea2t, sn2_15=sn2_15,
                     clobber_fibermap = clobber_fibermap, utah=utah,
                     termverbose=termverbose)
        logger = initializeLogger(config)
        writeVersionInfo(config, logger)

        #    Find correct MJD to start on
        initializeMJD(config, logger)
        #    Create poll workers and initialize file counts
        pollWorkers = createPollWorkers(config, logger)
        if config.catchup or config.redo: redo(pollWorkers, config, logger)
        else: runner(pollWorkers, config, logger)
        closeLogger(logger)

def parseNumList(string):
    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start,10), int(end,10)+1))
