#!/usr/bin/env python3
import boss_drp
from boss_drp.utils.splog import splog
from boss_drp.sos import sos_classes
from boss_drp.sos.sos_classes import SOS_config
from boss_drp.utils import putils, sxpar
from boss_drp.utils import sxpar
from boss_drp.utils.hash import create_hash
from boss_drp.prep.readfibermaps import readfibermaps
import boss_drp.sos.cleanup_sos  # This sets up cleanup for the main process
from boss_drp.sos.read_sos import read_SOS #log critical
splog._log.setLevel('DEBUG')
from boss_drp.sos.build_combined_html import build_combine_html
from boss_drp.sos import getSOSFileName
from boss_drp.sos.loadSN2Value import loadSN2Values

import functools
import builtins
import argparse
from argparse import ArgumentTypeError
import sys
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
import traceback

try:
    import sdssdb
    sdssdb.autoconnect = False
    from sdssdb.peewee.sdss5db.targetdb import database
    sdssdb.auto_reflect = False
    from sdssdb.peewee.sdss5db.targetdb import Design
    from sdssdb.peewee.sdss5db.opsdb import Configuration
    splog.add_external_handlers(sdssdb.log.name)
except:
    sdssdb = None

####
class fs_Config:
    """file sequencer Config Info"""

    def __init__(self, cfg):
        self.fitname  = ""
        self.fitdir   = ""
        self.plugname = ""
        self.plugdir  = ""
        self.plugging = ""
        self.flavor   = ""
        self.designMode = ""
        self.run_config = cfg

    def __str__(self):
        return ("fitname:    " + self.fitname + "\n" +
                "fitdir:     " + self.fitdir + "\n" +
                "plugname:   " + self.plugname + "\n" +
                "plugdir:    " + self.plugdir + "\n" +
                "plugging:   " + self.plugging + "\n" +
                "flavor:     " + self.flavor + "\n" +
                "designMode: " + self.designMode + "\n" +
                "------------------------------\n" +
                str(self.run_config));

def updateMJD(workers):
    """Check to see if a new MJD exists"""

    regex = sos_classes.Consts().MJDGlob;
    try:
        MJD = ls(SOS_config.fitsDir, regex)[-1][-5:]
        if (MJD == SOS_config.MJD):
            return

        SOS_config.MJD = MJD[-5:]
        for worker in workers:
            worker.fileCount = 0

        splog.info("Latest updated MJD found to be " + os.path.join(SOS_config.fitsDir, SOS_config.MJD))
    except:
        splog.critical("Could not find latest MJD in " + SOS_configs.fitsDir)
        splog.critical("GOODBYE!")
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


def Mode(cfg):
    start_time = time.time()
    if (sdssdb is not None) and (cfg.run_config.fps is True):
        try:
            if (not database.connected) or (not database.execute_sql("SELECT 1")):
                splog.info('connecting to sdssDB')
                try:
                    database.close()
                except:
                    pass
                database.connect()
            with database.atomic():
                dm = Design.select()\
                           .join(Configuration, on=(Configuration.design_id == Design.design_id))\
                           .where(Configuration.configuration_id == cfg.plugging)
            if len(dm) > 0:
                cfg.designMode = dm[0].design_mode_label
        except Exception as e:
            tb_str = traceback.format_exception(etype=type(e), value=e, tb=e.__traceback__)
            splog.critical("".join(tb_str))
        finally:
            pass
            #if database.connected:
            #    database.close()
    elif (not cfg.run_config.fps):
        cfg.designMode = 'Plate'
    else:
        cfg.designMode = 'Unknown'
    if cfg.designMode == '':
        cfg.designMode = 'Unknown'
    
    execution_time = time.time() - start_time
    splog.info(f"Execution time: {execution_time:.6f} seconds")
    return cfg
    
def previousExposure(cfg):
    """return a config for the previous exposure"""

    #    Parse: sdR-b1-00114186.fit.gz
    left  = cfg.fitname[:7]
    right = cfg.fitname[15:]
    exp   = (str(int(cfg.fitname[7:15]) - 1)).zfill(8)

    prevcfg = copy.copy(cfg)
    prevcfg.fitname = left + exp + right
   
    prvfitsExist = os.path.exists(os.path.join(prevcfg.fitdir, prevcfg.fitname))
    if prvfitsExist:
        prevcfg.plugging = plugging(prevcfg)
        prevcfg.flavor = flavor(prevcfg)
        prevcfg = Mode(prevcfg)
    splog.info("previous cfg:\n" + str(prevcfg))
    return prevcfg, prvfitsExist


def rule1(cfg):
    """Handle arc/flat ordering"""
    prevcfg, prvfitsExist = previousExposure(cfg)
    splog.info("Exposure Flavor: " + cfg.flavor)
    splog.info("Previous exposure exists: " + str(prvfitsExist))
    if prvfitsExist:
        splog.info("Previous Flavor: " + prevcfg.flavor)
        splog.info("Same Plugging: " + str(cfg.plugging == prevcfg.plugging))

    #    Handle flats -- process, and arc if was previous
    if cfg.flavor == "flat":
        splog.info("Exposure is a flat")
        processFile(cfg)
        if prvfitsExist and prevcfg.flavor  == "arc":
            splog.info("Processing previous arc")
            processFile(prevcfg)
        return True
    return False

def logecho(message, prefix=''):
    """Log a message and optionally print it based on termverbose."""
    # Add a prefix if provided
    message = f"{prefix}{message}" if prefix else message
    
    # Log the message
    splog.info(message)
    
    # Conditionally print the message if termverbose is True
    if SOS_config.termverbose:
        builtins.print(message)

class PrintRedirector:
    def __init__(self, logger_func):
        self.logger_func = logger_func
        self.original_print = builtins.print

    def __enter__(self):
        builtins.print = self.logger_func

    def __exit__(self, exc_type, exc_value, traceback):
        builtins.print = self.original_print

####
def processFile(cfg):
    """call sos_command on the file.  Will exit with error code if the command failts."""

    cmd  = "sos_command"
    cmd += " -f " + cfg.fitname
    cmd += " -i " + cfg.fitdir
    cmd += " -p " + cfg.plugname
    cmd += " -l " + cfg.plugdir
    cmd += " -s " + cfg.run_config.sosdir
    cmd += " -m " + cfg.run_config.MJD
    cmd += " -g " + cfg.designMode
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
    if cfg.run_config.bright_sn2:
        cmd += " -w "
    
    prefix = "sos_command(" + cfg.flavor + "): "
    logecho_wp = functools.partial(logecho, prefix=prefix)

    echo = cfg.run_config.termverbose

    i = 0
    license_crash = False
    
    from boss_drp.sos import filecheck
    ff = os.path.join(cfg.fitdir,cfg.fitname)

    while i < 5:
        if i > 0:
            time.sleep(2)
            license_crash = True
            splog.info("Trying again to get idl license")
        logecho("executing: " + cmd)
        if not filecheck.excellent(ff):
            ql = sxpar.sxparRetry(ff, "QUALITY", retries = 5)[0]
            logecho_wp(f'{ff} is not an excellent exposure ({ql})')
            continue
        else:
            logecho_wp(f'{ff} is an excellent exposure')
        if cfg.run_config.pause:
            logecho_wp(f"Sleeping for {sos_classes.Consts().licensePause} sec before "+
                       f"starting {cfg.fitname} reduction due to IDL License")
            time.sleep(sos_classes.Consts().licensePause)
        rv = putils.runCommand(cmd, echo=echo, prefix=prefix, logCmd=splog.info, limit=10)
        if rv[0] != 0:
            splog.info("\nCommand failed with rc = " + str(rv[0]) + "\n")
            sys.exit(1)
        if 'Failed to acquire license.' not in rv[1]:
            license_crash = False
            break
     
    if not license_crash:
        retry(postProcessFile, retries = 5, delay = 2, logger=splog.info, cfg=cfg)
        
    test = create_hash(os.path.join(cfg.run_config.sosdir,cfg.run_config.MJD))
    if test:
        splog.info("\nsha1sum is locked")

        
def postProcessFile(cfg):
    """call post sos_command commands on the file.  Will exit with error code if the command failts."""
    prefix = "sos_post_command(" + cfg.flavor + "): "
    logecho_wp = functools.partial(logecho, prefix=prefix)
    
    if cfg.flavor.lower() != 'science':
        logecho_wp(os.path.join(cfg.fitdir,cfg.fitname)+' is not a science frame')

        if cfg.flavor.lower() == 'arc':
            if (cfg.run_config.arc2trace) or (cfg.run_config.forcea2t):
                os.environ['BOSS_SPECTRO_REDUX'] = os.path.join(cfg.run_config.sosdir,f'{cfg.run_config.MJD}')
                
                cmd = (f"boss_arcs_to_traces --mjd {cfg.run_config.MJD} --no_hash "+
                       f"--obs {os.getenv('OBSERVATORY')} --cams {cfg.run_config.CCD} "+
                       f"--vers sos --threads 0 --sosdir {cfg.run_config.sosdir} "+
                       f"--fitsname {cfg.fitname}")
                logecho_wp(cmd)
                rv = putils.runCommand(cmd, echo=cfg.run_config.termverbose,
                                       prefix=prefix, logCmd=splog.info, env=os.environ.copy())
    else:
        sciE = getSOSFileName(os.path.join(cfg.fitdir,cfg.fitname))

        logecho_wp( os.path.join(cfg.fitdir,cfg.fitname)+' is a science frame')
        
        #load SN2 Values to DB
        with PrintRedirector(logecho_wp):
            logecho_wp( f'loadSN2Value -uv {os.path.join(cfg.run_config.sosdir,sciE)} {os.path.join(cfg.plugdir, cfg.plugname)}')
            loadSN2Value(os.path.join(cfg.run_config.sosdir,sciE),
                         os.path.join(cfg.plugdir, cfg.plugname),
                         verbose=True, update = True, sdssv_sn2=False)
            # read SOS
            logecho_wp( f'read_sos {cfg.run_config.sosdir} {cfg.run_config.MJD} --no_hash --exp={sciE}')
            warnings.warn = splog._original_warn # surpress the warning capture (read_SOS produces a lot of hidden warnings that I don't want to capture)
            read_SOS(cfg.run_config.sosdir, cfg.run_config.MJD, exp=sciE)
            warnings.warn = splog.Warning #turn warning capture back on
     
    with PrintRedirector(logecho_wp):
        # Build Index
        if cfg.run_config.CCD in ['b1','b2']:
            logecho_wp( f'build_combined_html {cfg.run_config.sosdir}')
            build_combine_html(cfg.run_config.sosdir, force=False)
    return



def sos_filesequencer(fitname, fitpath, plugname, plugpath):
    """
        Checks if processing a flat if there was an arc before the flat,
            then process the after the flat
        Otherwise process the fits file
    """


    config = fs_Config(SOS_config)
    config.fitname  = fitname
    config.fitdir   = fitpath
    config.plugname = plugname
    config.plugdir  = plugpath
    config.plugging = plugging(config)
    config.flavor = flavor(config)
    config = Mode(config)
    #    Rules return true if they processed the file and processing should stop
    #    process rule

    splog.info("Checking Rule 1")
    if not rule1(config):
        splog.info("Passed Rules; let's go!")
        processFile(config)

####
def processNewBOSSFiles(worker, files):
    """  Process new fits files

    Get the plugmap name and then add the appropiate sos command to the
    correctly numbered process list.

    Before the files are processed, they are sorted by name.  We really want the files
    sorted by time, but because of the sequence number, name is the same as time for any
    given camera.

    """

    #   Sort files by name to get into the right time order
    files.sort()
    splog.info("Sorted file list" + str(files))
    if len(files) == 0:
        splog.info('No Files')
    for file in files:
        splog.info("\n+++++++++++++++++++++++++++++++++++++++++++ \nprocessing new file: " + file)

        #- Get platetype from header (missing=BOSS)
        try:
            hdr = getheader(file)
        except:
            time.sleep(5)
            hdr = getheader(file)
            

        if 'FLAVOR' not in hdr:
            splog.info("Skipping exposure with missing FLAVOR keyword.")
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
            plugpath = getPlugMap(file)

            SOS_opts = {'confSummary':plugpath,'ccd':os.path.basename(file).split('-')[1], 'mjd':SOS_config.MJD,
                        'log':True, 'log_dir': os.path.join(SOS_config.sosdir,str(SOS_config.MJD))}
            readfibermaps(topdir=os.path.join(SOS_config.sosdir,str(SOS_config.MJD)), SOS=True, SOS_opts=SOS_opts,
                          logger=splog, clobber=SOS_config.clobber_fibermap)

            qf = os.path.abspath(file)
            fitname  = os.path.basename(qf)
            fitpath = os.path.dirname(qf)

            plugPath = os.path.abspath(plugpath)
            plugname  = os.path.basename(plugPath)
            plugpath = os.path.dirname(plugPath)

            sos_filesequencer(fitname, fitpath, plugname, plugpath)
        else:
            #- Don't crash if hdr is mangled and doesn't have EXPOSURE
            if 'EXPOSURE' in hdr:
                splog.info("Skipping %s exposure %d." % (platetype, hdr['EXPOSURE']))
            else:
                splog.info("Skipping %s exposure." % platetype)
            return

def writeVersionInfo():
    """Write a version string to a file"""

    verFile = os.path.join(SOS_config.controlDir, sos_classes.Consts().versionFile)
    f = open(verFile, "w")
    f.write(time.ctime() + " " + boss_drp.__version__ + "\n")
    f.close()
    splog.info("Version is %s" % boss_drp.__version__)
    splog.info("IDLSPEC2D Module is %s" % os.getenv('IDLSPEC2D_VER'))

####
def getPlugMap(file):
    """
    Returns the fully qualified name of the plugmap file
    """

    speclogDir = SOS_config.plugDir


    try:
        obs = os.environ['OBSERVATORY'].lower()
    except:
        splog.critical('Not running at APO or LCO, set OBS environmental variable to run')
        return ""
    plugmapDir = os.path.join(speclogDir, obs, 'summary_files')

    #   Get plugmap used by file
    try:
        plugmapFullId = sxpar.sxparRetry(file, "CONFID", retries = 5)[0]
    except TypeError as t:
        splog.critical("\nCould not parse " + file + "\n ->" + str(t))
        return ""
    except IndexError:
        splog.critical("\nKeyword CONFID not found in " + file)
        return ""
    try:
        if int(plugmapFullId) == -999:
            plugmapFullId = '0'
    except:
        plugmapFullId = '0'

    try:
        configgrp = '{:0>3d}XXX'.format(int(plugmapFullId)//1000)
        configdir = '{:0>4d}XX'.format(int(plugmapFullId)//100)
        plugmapDir = os.path.join(plugmapDir, configgrp, configdir)
    except:
        configgrp = '{:0>3d}XXX'.format(int(0)//1000)
        configdir = '{:0>4d}XX'.format(int(0)//100)
        plugmapDir = os.path.join(plugmapDir, configgrp, configdir)

    splog.info("Current confSummary directory is " +  plugmapDir)

    #   Parse plugmap name
    plugmapName   = "confSummaryF-" + plugmapFullId + ".par"
    plugParse     = plugmapFullId.split("-")
    plugmapId     = plugmapFullId #plugParse[1]
    splog.debug(file + " uses confSummaryF " + plugmapFullId + " with Id " + plugmapId)
    splog.debug("  full name of confSummaryF file is " + plugmapName)
    splog.debug("confId=" + plugmapId)# + ", pMJD=" + plugmapMJD)

    #   Check if the file exists, if not get it and add it to svn
    plugpath = os.path.join(plugmapDir, plugmapName)
    if os.path.isfile(plugpath):
        splog.info("Found existing confSummaryF file: " + plugpath)
    else:
        plugmapName   = "confSummary-" + plugmapFullId + ".par"
        plugParse     = plugmapFullId.split("-")
        plugmapId     = plugmapFullId #plugParse[1]
        splog.info("No confSummaryF file: " + plugpath)
        splog.debug(file + " uses confSummary " + plugmapFullId + " with Id " + plugmapId)
        splog.debug("  full name of confSummary file is " + plugmapName)
        splog.debug("confId=" + plugmapId)# + ", pMJD=" + plugmapMJD)

        plugpath = os.path.join(plugmapDir, plugmapName)

        if os.path.isfile(plugpath):
            splog.info("Found existing confSummary file: " + plugpath)
        else: splog.critical("Could not get confSummary for Id " + plugmapId)
    return os.path.abspath(plugpath)

def initializePollWorkers(workers):
    """Initialize poll workers with latest file counts"""

    for worker in workers:
        worker.fileCount = len(glob.glob(os.path.join(SOS_config.fitsDir, SOS_config.MJD, worker.glob)))
        splog.debug("\nInitialized PollWorker:\n" +  str(worker))

def initializeMJD():
    """Find the correct MJD to start looking for new files.  If the user specifies an MJD just test
    to see if it exists, otherwise, use the latest MJD."""

    #   First check for user specified
    if SOS_config.MJD != "0":
        path = os.path.join(SOS_config.fitsDir, SOS_config.MJD)
        if not os.path.isdir(path):
            splog.critical("Could not find user specified MJD path: " + path)
            splog.critical("GOODBYE!")
        splog.info("Using user specified MJD " + path)
    else:
        regex = sos_classes.Consts().MJDGlob;
        try:
            splog.debug("Looking for initial MJD in " + SOS_config.fitsDir)
            SOS_config.MJD = ls(SOS_config.fitsDir, regex)[-1][-5:]
            splog.info("Latest initial MJD found to be " + os.path.join(SOS_config.fitsDir, SOS_config.MJD))
        except:
            splog.critical("Could not find the latest MJD path: " + SOS_config.fitsDir)
            splog.critical("GOODBYE!")

def createPollWorkers():
    """Create poll workers"""

    workers = []

    num = 1
    for glob in SOS_config.globs:
        p = sos_classes.PollWorker()
        p.glob = glob
        p.workerNumber = num
        num += 1
        workers.append(p)
        splog.debug("\nnew PollWorker:\n" + str(p))

    return workers


def initializeLogger():
    """Startup logger and set the level"""

    lname = os.path.join(SOS_config.logDir, sos_classes.Consts().logName)
    if SOS_config.iname != "":
        lname += "-" + SOS_config.iname

    splog.set_SOS(sos_classes.Consts().logName,lname,SOS_config)
    splog.open()

    splog.info("Hello. " + sys.argv[0] + " started.")
    splog.info("Startup Configuration is: \n\n" + str(SOS_config) + "\n\n")

    return

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


def redo(workers):
    """Redo the command for files in the specified MJD"""

    #    Get files
    for worker in workers:
        files = lsltr(os.path.join(SOS_config.fitsDir, SOS_config.MJD), worker.glob)
        if SOS_config.exposure != None:
            allfiles = files
            files = []
            #    should only be one, but I do it this way to be sure things are working.
            #   also, I do the numeric check to avoid worrying about leading zeroes.
            for file in allfiles:
                splog.info("Checking exposure number of: " + file)
                exp = re.search("sdR\-..-(\d{8})\.fit.*$", file)
                if exp != None:
                    if int(exp.group(1)) == int(SOS_config.exposure):
                        splog.info("correct exposure number")
                        files.append(file)
        new = len(files)
        splog.info("Found " + str(new) + " files in " +
                  os.path.join(SOS_config.fitsDir, SOS_config.MJD, worker.glob))
        processNewBOSSFiles(worker, files)


def runner(pollWorkers):
    """
    Run the check for fits files and processes it
    """
        #   Create poll workers and initialize file counts
    initializePollWorkers(pollWorkers)
    #   Watch for new files.  Forever...  Unless there are exceptions.  Then
    #   try up to 5 times to get it working.  But mostly...  Forever!
    crashes = 5
    while crashes > 0:
        try:
            watch(pollWorkers)
        except SystemExit:
            raise
        except:
            crashes = crashes - 1
            if crashes > 0:
                splog.exception("!!! Uncaught exception in watch() !!!  Will Retry  !!!")
                time.sleep(5)
            else:
                splog.exception("!!! TOO MANY Uncaught exceptions in watch() !!!")
                raise

def watch(workers):
    """  Watch for new files

    When a new file comes in read the header to look for the plugmap and then check to see
    if the plugmap file exists.  If it doesn't, get the plugmap from the database and put it
    into the proper MJD directory.  Create the proper MJD directory for the plugmap if needed.

    Next, check to see if a newer MJD has been created.  If there are no new files and no new
    MJD then sleep for SOS_config.pollDelay.

    Note that only the latest MJD is ever checked, so once a new MJD is created only that MJD
    will be checked.
    """
    tpause = 0
    vpollDelay = 300
    disconnect = 3600
    spause = 0
    while True:
        pause = True
        #   First check for new files
        for worker in workers:
            files = lsltr(os.path.join(SOS_config.fitsDir, SOS_config.MJD), worker.glob)
            if len(files) != worker.fileCount:
                pause = False
                new = len(files) - worker.fileCount
                splog.info("Found " + str(new) + " new files in " +
                          os.path.join(SOS_config.fitsDir, SOS_config.MJD, worker.glob))
                #   File could get deleted...
                if new > 0:
                    processNewBOSSFiles(worker, files[-1 * new:])
                worker.fileCount = len(files)
                tpause = 0
                spause = 0
                disconnect = 3600

        #   Next check for a new MJD.  Don't wait if there's a new MJD
        if updateMJD(workers):  pause = False

        #   Pause if asked
        if pause:
            if (tpause % vpollDelay == 0) or (tpause / vpollDelay >= 1):
                splog.info(f"Sleeping for {SOS_config.pollDelay} seconds. (outout every {vpollDelay}s)")
                tpause = 0
            if (spause >= disconnect) and (disconnect != -1):
                try:
                    database.close()
                    disconnect = -1
                except:
                    pass
            spause += SOS_config.pollDelay
            tpause += SOS_config.pollDelay
            time.sleep(SOS_config.pollDelay)

def SOS(CCD, exp=None, mjd=None, catchup=False, redoMode=False,systemd=False, nodb=False,
        no_gz=False, no_reject=False, clobber_fibermap=False, sdssv_sn2=False,
        arc2trace=False, forcea2t=False, pause = False, test=False, utah=False,
        termverbose = False, sn2_15=False, bright_sn2=False, unlock = False):
    """
    The SOS controller for both manual runs and systemd tasks
    """
    try:
        for i, ex in enumerate(exp):
            if i > 1:
                pause = False
            #config = #sos_classes.Config();
            SOS_config.setup(CCD = CCD, mjd=mjd, exp=ex,
                             redo=redoMode, catchup=catchup, test=test, systemd=systemd,
                             no_gz=no_gz, nodb=nodb, no_reject=no_reject, sdssv_sn2=sdssv_sn2,
                             pause=pause, arc2trace=arc2trace, forcea2t=forcea2t, sn2_15=sn2_15,
                             clobber_fibermap = clobber_fibermap, utah=utah, bright_sn2=bright_sn2,
                             termverbose=termverbose)
            boss_drp.sos.cleanup_sos.check(force_unlock=unlock)
            initializeLogger()
            writeVersionInfo()

            #    Find correct MJD to start on
            initializeMJD()
            #    Create poll workers and initialize file counts
            pollWorkers = createPollWorkers()
            if SOS_config.catchup or SOS_config.redo: redo(pollWorkers)
            else: runner(pollWorkers)
            splog.close()
    except KeyboardInterrupt:
        splog.warning(f"SOS for CCD {CCD} interrupted in process {os.getpid()}.")
    except Exception as e:
        splog.warning(f"An error occurred in SOS for CCD {CCD}: {e}")
    finally:
        boss_drp.sos.cleanup_sos.cleanup()
        splog.close()
        
def parseNumList(string):
    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start,10), int(end,10)+1))
