#!/usr/bin/env python3

from astropy.io import fits
import os
import sys
import time
import getopt
from astropy.time import Time
from sdssdb.peewee.sdss5db import opsdb, targetdb
import sdssdb
print('sdssdb version:'+sdssdb.__version__)

####
class Config:
    """Config Info"""
    
    def __init__(self):
        self.verbose = False
        self.update  = False
        self.fits    = None
        self.confSum = None
        self.confID  = None
        self.design  = None
        self.sdssv_sn2 = False
    def __str__(self):
        if self.sdssv_sn2 is True:
            return ("Verbose:     " + str(self.verbose) + "\n" +
                    "Update:      " + str(self.update) + "\n" +
                    "fits:        " + self.fits + "\n" +
                    "confSummary: " + self.confSum + "\n" +
                    "confid:      " + self.confID+ "\n" +
                    "design:      " + self.design+ "\n" +
                    "sdssv_sn2:   " + self.sdssv_sn2+ "\n");
        else:
            return ("Verbose:     " + str(self.verbose) + "\n" +
                    "Update:      " + str(self.update) + "\n" +
                    "fits:        " + self.fits + "\n" +
                    "confSummary: " + self.confSum + "\n" +
                    "confid:      " + self.confID+ "\n" +
                    "design:      " + self.design+ "\n");


####
def usage():
    """Display usage and exit"""
    
    usageCMD = os.path.basename(sys.argv[0])

    print("usage:")
    print("\t%s [--update] [-v] fits-file confSummary-file" % usageCMD)
    print(" ")
    print("The fits file is the science frame output from apo-reduce.")
    print(" ")
    print("An error will occur if the exposure has already been processed, unless")
    print("--update is specified.")
    print(" ")
    print("Add -v for verbose output")

    sys.exit(1)
    

####
def parseCmdLine(args):
    """Parse command line arguments and return a Config"""

    # parse with options
    try:
        opts, pargs = getopt.gnu_getopt(sys.argv[1:], "v", ["update", "sdssv_sn2"])
    except:
        usage()

    if len(pargs) != 2:
        print("wrong number of files.\n")
        usage()

    cfg         = Config()
    cfg.fits    = pargs[0]
    cfg.confSum = os.path.basename(pargs[1])
    cfg.confID  = getConfig(cfg)
    cfg.design  = getDesign(cfg)

    #   Fill in the config
    for (opt, value) in opts:
        if opt == "-v":
            cfg.verbose = True
        if opt == "--update":
            cfg.update = True
        if opt == "sdssv_sn2":
            cfg.sdssv_sn2 = True
    #   Display config values if verbose
    if (cfg.verbose):
        print("Config values: \n" + str(cfg))
        
    return cfg


####
def getSN2(cfg):
    """Return the sn2 of the fits file as a string.  The SN2 is stored in FRAMESN2"""
    
    try: sn2list = str(fits.getval(cfg.fits, 'FRAMESN2'))
    except:
        print("WARNING: " + cfg.fits + " does not contain the header keyword FRAMESN2.")
        return 0.0;
    if cfg.sdssv_sn2:
        try: sn2list_v2 = str(fits.getval(cfg.fits, 'FSN2_v2'))
        except:
            print("WARNING: " + cfg.fits + " does not contain the header keyword FSN2_v2.")
            return 0.0;
 
    if cfg.verbose:
        print("(s/n)^2 is " + sn2list)
        if cfg.sdssv_sn2 is True:
            print("(s/n)^2 is " + sn2list_v2)

    if cfg.sdssv_sn2 is False:
        return sn2list
    else:
        return [sn2list, sn2list_v2]
    

####
def getMJD(cfg):
    """Return the MJD of the fits file as a string.  The MJD is taken from the fits header"""

    try: mjdlist = str(fits.getval(cfg.fits, 'MJD'))
    except:
        print("WARNING: " + cfg.fits + " does not contain the header keyword MJD.")
        return "00000"

    if cfg.verbose:
        print("MJD is " + mjdlist)

    return mjdlist


####
def getStartTime(cfg):
    """Return the start time of the fits file as a string.  The start time is taken from the fits header"""

    try: start2list = str(fits.getval(cfg.fits, 'TAI-BEG'))
    except:
        print("WARNING: " + cfg.fits + " does not contain the header keyword TAI-BEG.")
        return "0000000000.0"

    if cfg.verbose:
        print("TAI-BEG is " + start2list)
    start2list=str(float(start2list)/(24.0*3600.0))
    return start2list
    

####
def getExposureTime(cfg):
    """Return the exposure time of the fits file as a string.  The exposure time is taken from the fits header"""

    try: exptime2list = str(fits.getval(cfg.fits, 'EXPTIME'))
    except:
        print("WARNING: " + cfg.fits + " does not contain the header keyword EXPTIME.")
        return "0000.00"

    if cfg.verbose:
        print("EXPTIME is " + exptime2list)

    return exptime2list
    
def getConfig(cfg):
    try: cfid = str(fits.getval(cfg.fits, 'CONFID'))
    except:
        print("WARNING: " + cfg.fits + " does not contain the header keyword CONFID.")
        return("000000")
    if cfg.verbose:
        print("CONFID is " + cfid)
    return cfid

def getDesign(cfg):
    try: dsid = str(fits.getval(cfg.fits, 'DESIGNID'))
    except:
        print("WARNING: " + cfg.fits + " does not contain the header keyword DESIGNID.")
        return("00000")
    if cfg.verbose:
        print("EXPTIME is " + cfid)
    return dsid

####
def getFitsInfo(cfg):
    """Return the (exposure number, camera name) of the fits file (as strings)"""
    
    #   Parse science frame name.  e.g. sci-3690-b1-00105034.fits
    sciParse     = os.path.basename(cfg.fits).split("-")
    sciConf      = str(sciParse[1])
    sciCamera    = str(sciParse[2])
    sciExposure  = str(sciParse[3].split(".")[0])
    
    if cfg.verbose:
        print(("science exposure # " + sciExposure + " for configureID # " + sciConf +
              " for camera " + sciCamera))
        
    return (sciExposure, sciCamera)


####
def addOrUpdateExposure(cfg):
    """Add or update an exposure record with the (s/n)^2"""

    try:
        db_config =  opsdb.Configuration.get_or_create(configuration_id=cfg.confID, design=cfg.design)[0]
    except:
        db_config = opsdb.Configuration.create(design=cfg.design)
                            

    useTime = Time(float(getStartTime(cfg)),format="mjd").datetime #getMJD(cfg)

    
    #   In order to deal with the exposure we need some keys into other tables
    try:
        (expNum, camName) = getFitsInfo(cfg)
        camera = opsdb.Camera.get(label=camName)
        db_flavor = opsdb.ExposureFlavor.get(pk=1)  # science
#        survey  = session.query(Survey).filter_by(label="BHM").one()
    except:
        print("Could not retrieve camera, survey or instrument\n\n")
        raise

    #   Check to see if the exposure exists
   
    try:
        db_exp = opsdb.Exposure.get_or_none(configuration=db_config,
                                              #start_time=getStartTime(cfg),
#                                       survey = survey
                                              exposure_no = expNum)
                                #              start_time=useTime,
                                #              exposure_time = getExposureTime(cfg),
                                #             exposure_flavor=db_flavor)
    except:
        print("Problem querying for Exposure! \n\n")
        raise
    else:
        #   We don't do updates unless asked!
        if not cfg.update:
            print("Exposure record already exists!")
            print("Use '--update' to update the (s/n)^2 entry of an existing record\n\n")
            raise LookupError("Exposure record already exists.  Use --update to update")
        if cfg.verbose:
            if db_exp is None: print("Found an existing exposure record")
    if db_exp is None:
        if cfg.verbose:
            print("Creating new Exposure record.")
        db_exp = opsdb.Exposure.create(configuration=db_config,
                                       start_time=useTime,
#                                       survey = survey
                                       exposure_no = expNum,
                                       comment = "",
                                       exposure_time = getExposureTime(cfg),
                                       exposure_flavor=db_flavor)
                                       
#   We always update the quality (it can change)
#    exposure.exposure_status_pk = getQualityPK(cfg)
#    depreciated with FPS

    #   Check to see if the camera-frame exists
    cframe = None
    try:
        cframe = opsdb.CameraFrame.get_or_none(exposure=db_exp,
                                               camera=camera)

    except:
        print("Problem querying for CameraFrame! \n\n")
        raise
    else:
        #   We don't do updates unless asked!
        if not cfg.update:
            print("CameraFrame record already exists!")
            print("Use '--update' to update the (s/n)^2 entry of an existing record\n\n")
            raise LookupError("CameraFrame record already exists.  Use --update to update")
        if cfg.verbose:
            if cframe is not None: print("Found an existing CameraFrame record")

    #   Create an CameraFrame if needed

    sn2val=getSN2(cfg)
    if cfg.sdssv_sn2 is True:
        sn2val_2 = sn2val[1]
        sn2val = sn2val[0]
        comment = f"{int(time.time())}({sn2val},{sn2val_2})"
    else:
        comment = f"{int(time.time())}({sn2val})"
        #comment = str(int(time.time()))+"("+str(sn2val)+")"
    if cframe is None:
        if cfg.verbose:
            print("Creating new CameraFrame record.")
#        cframe = opsdb.CameraFrame.create(exposure=db_exp, camera=camera,
#                                          sn2=sn2val,sn2_v2=sn2val_2, comment = comment)
        cframe = opsdb.CameraFrame.create(exposure=db_exp, camera=camera,
                                          sn2=sn2val, comment = comment)
    else:
        cframe.comment = comment + "," + cframe.comment
#        cframe.sn2_v2 = sn2val_2
        cframe.sn2 = sn2val
        cframe.save()



####
def main(args):
    """Handle database connection and call the main script"""
    
    session = None
    config  = None
    
    try:
        config  = parseCmdLine(sys.argv[1:])
        if config.confID == '000000': return
        if config.design == '00000': return
        #(expNum, camName) = getFitsInfo(config)
        #camera = opsdb.Camera.get(label=camName)
        #print(camera)
        #db_flavor = opsdb.ExposureFlavor.get(pk=1)
        #print(db_flavor)
        addOrUpdateExposure(config)
    except:
        return

    
#### Start of script

if __name__=='__main__':
    main(sys.argv[1:])

