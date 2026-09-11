#!/usr/bin/env python3
from boss_drp.Config import config
from boss_drp.utils.splog import splog
from boss_drp.prep.spplan.tools import (build_exps, write_traceplan,
                                        obsTrace_mjdstart as obs_mjdstart)
from boss_drp.field import Fieldtype, Field
from boss_drp.utils import get_dirs, merge_ranges

from os import getenv
import os.path as ptt
from glob import glob
import time
from pydl.pydlutils.yanny import read_table_yanny


SDSSCOREVersion = getenv('SDSSCORE_VER', default= '')

        
def spplanTrace(flib=False, obs = None, mjd =None, sav_dir=None, 
                include_hartmann=False, exclude_arc = False, **extra_kwds):
    
    logfile = config.pipe['plan.trace.traceplan_logfile']
    if logfile is not None:
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    splog.info('spPlanTarget started at '+time.ctime())

    if obs is None:
        obs = config.pipe['fmjdselect.obs'].upper()
    if isinstance(obs, (list,tuple)):
        obs = obs[0]
    lco = True if obs.upper() == 'LCO' else False
    if lco:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
        OBS = 'LCO'
    else:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_N'
        OBS = 'APO'

    if mjd is None:
        mjd = config.pipe['fmjdselect.mjd']
    if not flib:
        mjdrange = config.pipe['fmjdselect.mjdrange']
        if mjdrange is None:
            mjdrange = [[obs_mjdstart[OBS], None]]
        else:
            if not isinstance(mjdrange[0], (list,tuple)):
                mjdrange = [mjdrange]
        for i, mjdr in enumerate(mjdrange):
            if mjd is not None:
                if isinstance(mjd, (list,tuple)):
                    if len(mjd) > 0:
                        if mjdrange[i][0] is None:
                            mjdrange[i][0] = int(mjd[0])
                        elif int(mjd[0]) < int(mjdrange[i][0]):
                            mjdrange[i][0] = int(mjd[0])
                else:
                    if mjdrange[i][0] is None:
                        mjdrange[i][0] = int(mjd)
                    elif int(mjd) < int(mjdrange[i][0]):
                        mjdrange[i][0] = int(mjd)
    else:
        mjdrange = [[mjd, None]]

    mjdrange = merge_ranges(mjdrange)
    #-------------
    # Determine the top-level of the output directory tree
    topdir = config.pipe['general.BOSS_SPECTRO_REDUX']
    splog.info('Setting TOPDIR='+topdir)
    
    run2d = config.pipe['general.RUN2D']
    splog.info('Setting RUN2D='+run2d)

    #----------
    # Read environment variable for BOSS_SPECTRO_DATA for finding raw data files.


    rawdata_dir = getenv(BOSS_SPECTRO_DATA)
    if rawdata_dir is None:
        splog.info('ERROR: Must set environment variable BOSS_SPECTRO_DATA')
        exit(1)
    
    speclog_dir = getenv('SPECLOG_DIR')
    if speclog_dir is None:
        splog.info('ERROR: Must set environment variable SPECLOG_DIR')
        exit()
    splog.info('Setting SPECLOG_DIR='+speclog_dir)
    
    #----------
    # Create a list of the MJD directories (as strings)
    mjd_plans = config.pipe['fmjdselect.trace_all_mjds']
    if mjd_plans:
        splog.info('Limiting MJDs to those with existing spPlan2d files')
        mjd_plans = []
        fc = Field(topdir, run2d, '*')
        plans2d_tmp = glob(ptt.join(fc.dir(), 'spPlan2d*'))
        for plan2d in plans2d_tmp:
            if ptt.basename(plan2d).split('.')[0].split('-')[-1] in mjd_plans:
                continue
                
            plan = read_table_yanny(plan2d,'SPEXP')
            if plan.meta['OBS'] == OBS:
                mjd_plans.append(str(plan.meta['MJD']))
        mjd_plans = list(set(mjd_plans))
        if mjd is not None:
            mjd = list(set(mjd_plans) & set(mjd))
        else:
            mjd = mjd_plans


    mjdlist = get_dirs(rawdata_dir, subdir='', pattern='*', match=mjd, ranges = mjdrange)
    nmjd = len(mjdlist)
    splog.info(f'Number of MJDs = {nmjd}')
    if nmjd == 0:
        splog.info('No Valid MJDs')
        return None
    plateflavors = ['BHM', 'BHM&MWM', 'EBOSS', 'BOSS']
    #---------------------------------------------------------------------------
    # Loop through each input MJD directory

    legacy = config.pipe['SDSS_Generation.legacy']
    plates = config.pipe['SDSS_Generation.plates']
    fps = config.pipe['SDSS_Generation.fps']
    if config.pipe['SDSS_Generation.sdssv']:
        plates = True
        fps = True
    no_remote = not config.pipe['general.REMOTE']
    release = config.pipe['general.RELEASE']
    verbose = config.pipe['plan.daily.traceplan_verbose']
    clobber = config.pipe['Clobber.clobber_spTrace']
    override_manual = config.pipe['pipe.trace.override_manual_trace']
    dithered_pmjds = []
    for i, mj in enumerate(mjdlist):
        ftype = Fieldtype(fieldid=None, mjd=mj, obs=OBS)
        if not ftype.check(legacy=legacy, plates=plates, fps=fps):
            continue
        splog.info('----------------------------')
        splog.info(f'MJD: {mj} {ftype} ({i+1} of {len(mjdlist)})')
        allexps, ftype = build_exps(i, mj, mjdlist, OBS, rawdata_dir, ftype, spplan_Trace=True,
                                    legacy=legacy, plates=plates, fps=fps, lco=lco, 
                                    include_hartmann=include_hartmann, exclude_arc = exclude_arc,
                                    no_remote=no_remote, release=release, verbose=verbose)
        write_traceplan(allexps, mj, ftype, OBS, exclude_arc, sav_dir, topdir, 
                        run2d, release, clobber, override_manual)

        del allexps
    splog.info('----------------------------')
    splog.info('Successful completion of spplanTrace at '+ time.ctime())

    if logfile is not None:
        splog.close()
    return(nmjd)
        
