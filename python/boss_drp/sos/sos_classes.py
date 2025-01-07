#!/usr/bin/env python

import os, sys
import os.path as ptt

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
        self.licensePause = "5"
    
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
        self.fps        = True
        self.exposure   = None
        self.fitsDir    = '/data/spectro'
        self.plugDir    = os.getenv('SDSSCORE_DIR')
        self.controlDir = '/home/sdss5/boss/sos/control'
        self.logDir     = '/home/sdss5/boss/sos/logs'
        self.sosdir     = '/data/boss/sos'
        self.dlogLevel  = 40
        self.logLevel   = 40
        self.iname      = ""
        self.globs      = ["*"]
        self.command    = ""
        self.pollDelay  = 2
        self.nosvn      = True
        self.nocal      = True
        self.bookkeep   = False
        self.nice       = False
        self.platedb    = False
        self.redo       = False
        self.catchup    = False
        self.utah       = False
        self.test       = False
        self.nodb       = False
        self.no_reject  = False
        self.clobber_fibermap = False
        self.sdssv_sn2  = False
        self.sn2_15     = False
        self.arc2trace  = False
        self.forcea2t   = False
        self.pause      = False
        self.verbose    = 6
        self.termverbose=False
        self.set_logLevel()
        
    def set_logLevel(self):
        #   Don't want to apply -v on each call, so always start with a base
        if (self.verbose > 0):
            self.logLevel = max(1, self.dlogLevel - self.verbose * 10)
            
    def setup(self,  CCD='b1', mjd = None, exp = None,
              redo=False, catchup=False, test=False, systemd=False,
              no_gz=False, nodb=False, no_reject = False, sdssv_sn2 = False,
              pause = False, arc2trace=False, forcea2t=False, sn2_15 = False,
              clobber_fibermap=False, utah=False, termverbose=False,
              bright_sn2 = False):
        self.nodb = nodb
        self.no_reject = no_reject
        self.sdssv_sn2 = sdssv_sn2
        self.sn2_15 = sn2_15
        self.arc2trace = arc2trace
        self.forcea2t = forcea2t
        self.pause = pause
        self.clobber_fibermap = clobber_fibermap
        self.set_logLevel()
        self.iname = CCD
        self.termverbose = termverbose
        self.bright_sn2 = bright_sn2
        if bright_sn2:
            self.sn2_15 = True

        if no_gz:
            self.globs=["sdR-"+CCD+"-*.fit.gz","sdR-"+CCD+"-*.fit"]
        else:
            self.globs=["sdR-"+CCD+"-*.fit.gz"]

        if systemd:
            self.redo = False
            self.catchup = False
            exp = None
            mjd = None
            self.termverbose = False
            
        if redo:
            self.sosdir='/data/boss/sosredo'
            self.nice=True
            self.redo=True
            self.iname = self.iname+'_redo'
        elif test:
            self.sosdir='/data/boss/sosredo/dev'
            self.nice=True
            self.redo=True
            self.iname = self.iname+'_dev'
            self.logDir = '/data/boss/sosredo/dev/logs/'
            self.controlDir = '/data/boss/sosredo/dev/control/'
            self.nodb = True
            self.test = True
        elif catchup:
            self.nice=True
            self.redo=True
            self.catchup = True
            self.iname = self.iname+'_catchup'
        elif utah:
            self.redo=True
            self.nodb=True
            self.utah=True
            self.iname = self.iname+'_utah'
            self.nice=True
            self.pause=False
            self.sosdir = ptt.join(os.getenv('BOSS_SPECTRO_REDUX'), 'sos', os.getenv('RUN2D', default=''), '{obs}')
            self.sosdir = self.sosdir.format(obs=os.getenv('OBSERVATORY').lower())
            self.logDir = ptt.join(self.sosdir, 'logs')
            self.controlDir = ptt.join(self.sosdir, 'control')
            if os.getenv('OBSERVATORY').lower() == 'apo':
                self.fitsDir = os.getenv('BOSS_SPECTRO_DATA_N')
            else:
                self.fitsDir = os.getenv('BOSS_SPECTRO_DATA_S')
        else:
            pass
            
        if exp is not None:
            self.exposure = str(exp)
            self.iname = self.iname+'_'+str(exp).zfill(8)
        if mjd is not None:
            self.MJD = mjd
            self.iname = self.iname+'_'+str(mjd)

        if not ptt.exists(self.logDir):
            os.makedirs(self.logDir, exist_ok = True)
        if not ptt.exists(self.controlDir):
            os.makedirs(self.controlDir, exist_ok = True)
        if not ptt.exists(self.sosdir):
            os.makedirs(self.sosdir, exist_ok = True)

            
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return ("Config:\n"    +
                "MJD:        " + self.MJD + "\n" +
                "Exposure    " + str(self.exposure) + "\n" +
                "sosDir      " + self.sosdir + "\n" +
                "fitsDir:    " + self.fitsDir + "\n" +
                "plugDir:    " + self.plugDir + "\n" +
                "controlDir: " + self.controlDir + "\n" +
                "logDir:     " + self.logDir + "\n" +
                "logLevel:   " + str(self.logLevel) + "\n" +
                "FPSMode:    " + str(self.fps) + "\n" +
                "iname:      " + self.iname + "\n" +
                "glob:       " + str(self.globs) + "\n" +
                "command:    " + self.command + "\n" +
                "pollDelay:  " + str(self.pollDelay) + "\n" +
                "nice:       " + str(self.nice) + "\n" +
                "plateDb:    " + str(self.platedb) + "\n" +
                "redo:       " + str(self.redo) + "\n" +
                "test:       " + str(self.test) + "\n" +
                "catchup:    " + str(self.catchup) + "\n" +
                "bookkeep:   " + str(self.bookkeep) + "\n" +
                "NoSvn:      " + str(self.nosvn) + "\n" +
                "noDB:       " + str(self.nodb) + "\n" +
                "no_reject:  " + str(self.no_reject) + "\n" +
                "clobber_fibermap: " + str(self.clobber_fibermap)+ "\n" +
                "sdssv_sn2:  " + str(self.sdssv_sn2) + "\n" +
                "sn2_15:     " + str(self.sn2_15) + "\n" +
                "arc2trace:  " + str(self.arc2trace) + "\n" +
                "forcea2t:   " + str(self.forcea2t) + "\n"+
                "pause:      " + str(self.pause));
                   

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
        
