#!/usr/bin/env python3
import boss_drp
from boss_drp.utils.splog import splog
from boss_drp.prep.spplan import (build_exps, get_master_cal, find_nearest, pair_ccds,
                                  obsTrace_mjdstart as obs_mjdstart, spplan_findrawdata,
                                  write_plan)
from boss_drp.field import field_to_string, Fieldtype, Field
from boss_drp.utils import (Sphdrfix, mjd_match, get_dirs, getcard)
from boss_drp.prep.GetconfSummary import find_confSummary, find_plPlugMapM, get_confSummary
from boss_drp.utils.reject import Reject

from sdss_access.path import Path
from sdss_access import Access
from sdss_access import __version__ as saver
from tree import __version__ as treever

from os import getenv, makedirs, rename
import os.path as ptt
from glob import glob
import time
from astropy.table import Table, vstack, Column, unique
from astropy.io import fits
from collections import OrderedDict
from pydl.pydlutils.yanny import read_table_yanny, yanny, write_table_yanny
from pydl import __version__ as pydlVersion
import subprocess
import numpy as np


SDSSCOREVersion = getenv('SDSSCORE_VER', default= '')
idlspec2dVersion = boss_drp.__version__

        
def spplanTrace(topdir=None, run2d=None, mjd=None, mjdstart=None, mjdend=None,
             lco=False, clobber=False, release='sdsswork', logfile=None, no_remote=True,
             legacy=False, plates=False, fps = True, override_manual=False, sav_dir=None,
             verbose = False, no_dither = False, mjd_plans=False, flib=False, **extra_kwds):
    
    if logfile is not None:
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    splog.info('spPlanTarget started at '+time.ctime())

    if lco:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
        OBS = 'LCO'
    else:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_N'
        OBS = 'APO'

    if not flib:
        if mjdstart is None:
            mjdstart = obs_mjdstart[OBS]
            splog.info(f'{OBS} TraceFlat not valid before {mjdstart}... Setting mjdstart = {mjdstart}')
        elif int(mjdstart) < obs_mjdstart[OBS]:
            mjdstart = obs_mjdstart[OBS]
            splog.info(f'{OBS} TraceFlat not valid before {mjdstart}... Resetting mjdstart = {mjdstart}')
    else:
        mjdstart = mjd
    mjdstart = int(mjdstart)
    #-------------
    # Determine the top-level of the output directory tree
    if topdir is None:
        topdir = getenv('BOSS_SPECTRO_REDUX')
    splog.info('Setting TOPDIR='+topdir)
    
    if run2d is None:
        run2d = getenv('RUN2D')
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
    if mjd_plans:
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


    mjdlist = get_dirs(rawdata_dir, subdir='', pattern='*', match=mjd, start=mjdstart, end=mjdend)
    nmjd = len(mjdlist)
    splog.info(f'Number of MJDs = {nmjd}')
    if nmjd == 0:
        splog.info('No Valid MJDs')
        return None
    plateflavors = ['BHM', 'BHM&MWM', 'EBOSS', 'BOSS']
    #---------------------------------------------------------------------------
    # Loop through each input MJD directory

    dithered_pmjds = []
    for i, mj in enumerate(mjdlist):
        ftype = Fieldtype(fieldid=None, mjd=mj)
        if not legacy:
            if ftype.legacy is True:
                return None
        if not plates:
            if ftype.plates is True:
                return None
        if not fps:
            if ftype.fps is True:
                return None
        splog.info('----------------------------')
        splog.info(f'MJD: {mj} {ftype} ({i+1} of {len(mjdlist)})')
        allexps, ftype = build_exps(i, mj, mjdlist, OBS, rawdata_dir, ftype, spplan_Trace=True,
                                    legacy=legacy, plates=plates, fps=fps, lco=lco,
                                    no_remote=no_remote, release=release, verbose=verbose)
        thismjd = int(mj)
        if len(allexps) == 0:
            splog.info('No Calibration Frames for '+mj)
        else:
            if ftype.legacy or ftype.plates:
                fieldmap_col = 'mapname'
            else:
                fieldmap_col = 'fieldid'

            
            manual = 'F'
            allexps = get_master_cal(allexps, obs=OBS, mjd=mj)
            if allexps is None:
                continue
            planfile = 'spPlanTrace-' + mj + '_'+OBS+'.par'
            if sav_dir is None:
                planfile = ptt.join(topdir, run2d, 'trace', str(thismjd), planfile)
            else:
                planfile = ptt.join(sav_dir, str(thismjd), planfile)

            meta=OrderedDict({
                            'MJD':              mj                       +"   # Modified Julian Date",
                            'OBS':              OBS                      +"   # Observatory",
                            'RUN2D':            run2d                    +"   # 2D reduction name",
                            'idlspec2dVersion': "'"+idlspec2dVersion+"'" +"   # idlspec2d Version when building plan",
                            #'idlutilsVersion':  "'"+idlutilsVersion+"'"  +"   # idlutils Version when building plan",
                            'pydlVersion':      "'"+pydlVersion+"'"      +"   # Version of pydl when building plan",
                            #'speclogVersion':   "'"+speclogVersion+"'"   +"   # speclog Version when building plan",
                            'SDSSCOREVersion':  "'"+SDSSCOREVersion+"'"  +"   # SDSSCORE Version when building plan",
                            'SDSS_access_Ver':  "'"+saver+"'"            +"   # sdss_access Version when building plan",
                            'sdss_tree_Ver':    "'"+treever+"'"          +"   # sdss-tree Version when building plan",
                            'SDSS_access_Release': "'"+release+"'"       +"   # SDSS-access Release Version when building plan",
                            'manual':            manual                  +"   # Manually edited plan file (T: True, F: False)"
                                    })
            allexps = pair_ccds(ftype, allexps, OBS=OBS)
            write_plan(planfile, allexps, meta=meta, clobber=clobber, override_manual=override_manual)
        del allexps
    splog.info('----------------------------')
    splog.info('Successful completion of spplanTrace at '+ time.ctime())

    if logfile is not None:
        splog.close()
    return(nmjd)
        

