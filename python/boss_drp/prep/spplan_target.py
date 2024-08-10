#!/usr/bin/env python3
import boss_drp
from boss_drp.utils import (match, load_env, Splog)
from boss_drp.field import field_dir
from boss_drp.utils import jdate

splog = Splog()
try:
    from sdssdb.peewee.sdss5db.targetdb import database
    test = database.set_profile(load_env('DATABASE_PROFILE', default='pipelines'))
    from sdssdb.peewee.sdss5db.targetdb import CartonToTarget, Carton, Target
    nodb=False
except:
    if load_env('DATABASE_PROFILE', default='pipelines').lower() in ['pipelines','operations']:
        splog.log('ERROR: No SDSSDB access')
        exit()
    else:
        splog.log('WARNING: No SDSSDB access')
    nodb=True
from sdss_semaphore.targeting import TargetingFlags

import numpy as np
from os import getenv, makedirs
import os.path as ptt
import sys
from astropy.table import Table, vstack
import astropy.time
import argparse
from collections import OrderedDict
from pydl.pydlutils import yanny
from subprocess import getoutput
from datetime import date
from glob import glob
from tqdm import tqdm
import time



idlspec2dVersion = boss_drp.__version__

plan_fill = {'TARGID':-1,      'FMJD_LIST':'', 'CATALOGID_LIST':-1,
             'FIELDS_LIST':-1, 'MJD_LIST':-1,  'EPOCH_COMBINE': -1}
plan_dtyp = {'TARGID':int,      'FMJD_LIST':object, 'CATALOGID_LIST': int,
            'FIELDS_LIST':int,  'MJD_LIST':int ,    'EPOCH_COMBINE': int}


def clean_nones(array):
    if array is not None:
        if isinstance(array, (list, tuple, np.ndarray)):
            array = np.asarray(array).tolist()
            while None in array:
                array.remove(None)
        if isinstance(array,(float,int,str)):
            array= np.asarray(array)
        if len(array) == 0:
            array = None
        else:
            array = np.asarray(array)
    return(array)

def read_SPALL(topdir, run2d, cartons=None, catalogids=None, program=None,
                use_catid=False, use_firstcarton=False, usedb=False):
    """
    Read in spall file and filter to only preserve targets matching criteria
    """
    splog.log('Reading spAll')
    if ptt.exists(ptt.join(topdir, run2d, 'summary','daily','spAll-'+run2d+'.fits.gz')):
        spAll = Table.read(ptt.join(topdir, run2d, 'summary','daily', 'spAll-'+run2d+'.fits.gz'), 1)
    elif ptt.exists(ptt.join(topdir, run2d, 'spAll-'+run2d+'.fits.gz')):
        spAll = Table.read(ptt.join(topdir, run2d, 'spAll-'+run2d+'.fits.gz'), 1)
    elif ptt.exists(ptt.join(topdir, run2d, 'spAll-'+run2d+'.fits')):
        spAll = Table.read(ptt.join(topdir, run2d, 'spAll-'+run2d+'.fits'), 1)
    else:
        splog.log(ptt.join(topdir, run2d, 'summary','daily','spAll-'+run2d+'.fits')+' is missing')
        exit()
    sel = spAll.copy()
    selected = Table()
    
    program = clean_nones(program)
    cartons = clean_nones(cartons)
    catalogids = clean_nones(catalogids)


    if program is not None:
        selected = Table()
        for prog in program:
            selected = vstack([selected,sel[np.where(match(np.char.decode(aspAll['PROGRAMNAME'].data), prog))[0]]])
        sel = selected.copy()
    if cartons is not None:
        if use_firstcarton:
            selected = Table()
            for carton in cartons:
                selected = vstack([selected,sel[np.where(match(np.char.decode(spAll['FIRSTCARTON'].data), carton))[0]]])
            sel = selected.copy()
        else:
            if nodb:
                splog.info('Reverting to useDB = False')
                usedb = False
            if usedb:
                selected = Table()
                catids = []
                for carton in tqdm(cartons,desc='Carton',position=1, leave=False, disable=True):
                    carton = carton.replace('*','')
                    tp = CartonToTarget.select(Target.catalogid).join(Carton).switch().join(Target)\
                                       .where(Carton.carton.contains(carton)).dicts() #.limit(1000)
                    for t in tqdm(tp, desc='Carton2targ', position=2, leave=False,disable=True):
                        catids.append(t['catalogid'])
                catids = list(set(catids))
                mask = np.in1d(sel['CATALOGID'].data, catids)
                selected = sel[np.where(mask)]
                sel = selected.copy()
            else:
                check_cartons = []
                flags = TargetingFlags(sel['SDSS5_TARGET_FLAGS'].data)
                for carton in cartons:
                    idx = np.where(match(flags.all_carton_names,carton))[0]
                    check_cartons.extend(np.asarray(flags.all_carton_names)[idx].tolist())
                masks = None
                for carton in check_cartons:
                    tm = flags.in_carton_name(carton)
                    if masks is None:
                        masks = tm
                    else:
                        masks = [any(tup) for tup in zip(masks, tm)]
                sel = sel[masks]
                selected = sel
    if catalogids is not None:
        selected = Table()
        cid_col = 'CATALOGID' if use_catid is True else 'SDSS_ID'
        for catid in catalogids:
            selected = vstack([selected, sel[np.where(sel[cid_col] == int(catid))[0]]])
    if len(selected) > 0:
        selected['FMJD']=np.char.add(np.char.add(np.char.zfill(selected['FIELD'].data.astype(str),6), '-'),selected['MJD'].astype(str))
    return(selected)


def mjd_filter(spAll, mjd=None, mjdstart=None, mjdend=None):
    """
    Filters spAll time to the mjds and/or mjd range probvided
    """
    selected=Table()
    mjd = clean_nones(mjd)
    mjdstart = clean_nones(mjdstart)
    mjdend = clean_nones(mjdend)
    
    if mjd is not None or mjdstart is not None or mjdend is not None:
        splog.log('Filtering MJDs')

    if mjd is not None:
        for mj in mjd:
            selected = vstack([selected, spAll[np.where(spAll['MJD'] == int(mj))[0]]])
        spAll=selected.copy()
    if mjdstart is not None:
        spAll=spAll[np.where(spAll['MJD'] >= int(mjdstart))[0]]
    if mjdend is not None:
        spAll=spAll[np.where(spAll['MJD'] <= int(mjdend))[0]]
    return(spAll)



def build_plan(spAll, use_catid=False, coadd_mjdstart = None):
    """
    Buld plan file from filtered spAll file
    """
    plan = Table()
    splog.log('Building spplan_target')
    cid_col = 'CATALOGID' if use_catid is True else 'SDSS_ID'
    catalogids = np.asarray(spAll[cid_col])
    catalogids = catalogids[np.where(catalogids != -999)[0]]
    catalogids = np.unique(catalogids)
    catid_cols = (np.asarray(spAll.colnames))[np.where(match(spAll.colnames, 'CATALOGID*'))[0]]
    for catid in tqdm(catalogids, desc=f'Building {cid_col} Plan', position=2, disable=False):
        idx    = np.where(spAll[cid_col] == catid)[0]
        fmjds  = np.asarray(spAll[idx]['FMJD'].data).astype(object)
        mjds   = spAll[idx]['MJD'].data
        fields = spAll[idx]['FIELD'].data
        epoch  = np.max(mjds)
        catids = []
        for cc in catid_cols:
            catids.extend((spAll[idx][cc].data).tolist())
        if coadd_mjdstart is not None:
            if epoch < int(coadd_mjdstart):
                splog.log(f'Skipping {catid} (EPOCH_COMBINE={epoch} < {coadd_mjdstart})')
                continue
        catids = list(set(catids))
        catids = [i for i in catids if i is not None]

        new = Table({'TARGID':[catid],'FMJD_LIST':[fmjds],
                     'FIELDS_LIST':[fields], 'MJD_LIST':[mjds],
                     'CATALOGID_LIST':[catids],'EPOCH_COMBINE':[epoch]})

        if len(plan) > 0:
            for col in new.colnames:
                if len(new[col].shape) == 1:
                    continue
                if new[col].shape[1] > plan[col].shape[1]:
                    pad = np.full([plan[col].shape[0], new[col].shape[1] - plan[col].shape[1]], plan_fill[col], dtype=plan_dtyp[col])
                    plan[col] = np.hstack([plan[col], pad])
                if new[col].shape[1] < plan[col].shape[1]:
                    pad = np.full([new[col].shape[0], plan[col].shape[1] - new[col].shape[1]], plan_fill[col], dtype=plan_dtyp[col])
                    new[col] = np.hstack([new[col], pad])
        
        plan = vstack([plan,new])
        
        #if len(plan)> 400: break
    if len(plan) == 0:
        splog.log('No Valid TARGETS')
        plan = None
    else:
        plan['FMJD_LIST']=plan['FMJD_LIST'].astype(str)
        plan = mod_epoch(plan)
    return(plan)

def mod_epoch(plan):
    epochs = np.sort(np.unique(plan['EPOCH_COMBINE'].data))
    EPOCH_COMBINE = plan['EPOCH_COMBINE']
    for i, ec in enumerate(epochs):
        idx = np.where(plan['EPOCH_COMBINE'].data == ec)[0]
        if len(idx) <= 10:
            if i + 1 == len(epochs):
                splog.log(f'Dropping EPOCH_COMBINE={ec}')
                EPOCH_COMBINE[idx] = 0
            else:
                EPOCH_COMBINE[idx] = epochs[i+1]
    plan['EPOCH_COMBINE'] = EPOCH_COMBINE
    idx = np.where(plan['EPOCH_COMBINE'].data != 0)[0]
    plan = plan[idx]
    return(plan)

def write_plan(name, plan, topdir, run2d, run1d, mjd_start, mjd_end, coadd_mjdstart = None,
               clobber=False, rerun1d=False, use_catid=False, obs=None):
    """
    Write plan to file
    """
    cid_col = 'CATALOGID' if use_catid is True else 'SDSS_ID'
    topdir2d = ptt.join(topdir, run2d)
    planfile = ptt.join(field_dir(topdir2d, name, custom=True),
                        'spPlanCustom-'+name+'-'+jdate.astype(str)+'.par')
    makedirs(ptt.dirname(planfile), exist_ok=True)
    plan.meta=OrderedDict({'NAME':             name                   +"  # Name of Custom Coadd",
                           'RUN2D':            run2d                  +"  # Run2d Version",
                           'RUN1D':            run1d                  +"  # Run1d Version",
                           'idlspec2dVersion': idlspec2dVersion       +"  # Version of idlspec2d when building plan file",
                           #'idlutilsVersion':  idlutilsVersion        +"  # Version of idlutils when building plan file",
                           'Rerun_RUN1D':      str(rerun1d)           +"  # 1D anylsis of Custom Coadd",
                           'Date':time.strftime("%m/%d/%Y-%H:%M",time.localtime())+"  # Date of creation",
                           'CreateMJD' :       jdate.astype(str)      +"  # MJD of creation",
                           'MJD_range':  '-'.join([str(mjd_start),str(mjd_end)])    +"  # Range of MJDs available",
                           'TARGID':           cid_col                +"  # TARGID column maps to "+cid_col,
                           'MJD':              jdate.astype(str)      +"  # MJD of Coadd"
                           })
    if obs is not None:
        plan.meta['OBS'] = obs.upper()
    if coadd_mjdstart is not None:
        plan.meta['MJDSTART'] = coadd_mjdstart
    if ptt.exists(planfile):
        if clobber is False:
            splog.info('WARNING: Will not over-write plan file: ' + ptt.basename(planfile))
            return
        else:
            splog.info('WARNING: Over-writing plan file: ' + ptt.basename(planfile))
    else:
        splog.info('Writing plan file '+ ptt.basename(planfile))
    plan.convert_unicode_to_bytestring()
    yanny.write_table_yanny(plan, planfile,tablename='COADDPLAN', overwrite=clobber)
    return

def CustomCoadd(name, topdir, run2d, run1d, cartons=None, catalogids=None, obs=None,
                clobber=False, logfile=False, mjd=None, mjdstart=None, mjdend=None,
                program=None, rerun1d=False, use_catid=False, use_firstcarton=False,
                coadd_mjdstart = None,usedb=False):
    """
    Run a CustomCoadd Schema
    """
    splog.open(logfile=logfile)
    if logfile is not None: splog.log('Log File '+ptt.basename(logfile)+' opened '+time.ctime())

    if topdir is None: topdir = getenv('BOSS_SPECTRO_REDUX')
    if run2d is None: run2d = getenv('RUN2D')
    if run1d is None: run1d = getenv('RUN1D')

    spAll = read_SPALL(topdir, run2d, cartons=cartons, catalogids=catalogids,
                        program=program, use_catid=use_catid,
                        use_firstcarton = use_firstcarton,
                        usedb = usedb)
    spAll = mjd_filter(spAll, mjd=mjd, mjdstart=mjdstart, mjdend=mjdend)
    spAll = zwarn_filter(spAll)
    if obs is not None:
        spAll = obs_filter(spAll, obs)
    if len(spAll) == 0:
        splog.log(f'No Targets matching {name}')
        return
    
    if mjdstart is None:
        mjdstart = np.min(spAll['MJD'].data)
    if mjdend is None:
        mjdend = np.max(spAll['MJD'].data)
    if mjd is not None:
        mjd_start = mjd_end = mjd
    plan = build_plan(spAll, use_catid=use_catid, coadd_mjdstart = coadd_mjdstart)
    if plan is not None:
        write_plan(name, plan, topdir, run2d, run1d, mjdstart, mjdend,clobber=clobber,
                    rerun1d=rerun1d, use_catid=use_catid, obs=obs,
                    coadd_mjdstart=coadd_mjdstart)
                
    splog.close()
    return

def obs_filter(spAll, obs):
    splog.log(f'Filtering to {obs} observations')
    idx = np.where(np.char.upper(spAll['OBS'].data.astype(str)) == obs.upper())[0]
    return(spAll[idx])

def zwarn_filter(spAll):
    splog.log('Filtering with ZWARNING (UNPLUGGED,BAD_TARGET,NODATA)')
    idx = np.where(spAll['ZWARNING'].data < 128)[0]
    return(spAll[idx])


def get_key(fp):
    filename = ptt.splitext(ptt.basename(fp))[0]
    return(int(filename.split('-')[-1]))



def batch(topdir, run2d, run1d, DR = False, clobber=False,logfile=None, coaddfile=None,
        name=None, use_catid=False, use_firstcarton=False, coadd_mjdstart = None, obs=None,
        usedb=False):
    """
    Batch Build all Schema provided by the coaddfile
    """

    if topdir is None: topdir = getenv('BOSS_SPECTRO_REDUX')
    if run2d is None: run2d = getenv('RUN2D')
    if run1d is None: run1d = getenv('RUN1D')
    if coaddfile is None:
        coaddfile = ptt.join(topdir,run2d,'fields','SDSSV_BHM_COADDS.par')

    COADDS = yanny.read_table_yanny(coaddfile, 'SCHEMA')
    COADDS.convert_bytestring_to_unicode()

    COADDS = COADDS[np.where(COADDS['ACTIVE'] == 1)[0]]
    
    COADDS['MJD'] = COADDS['MJD'].astype(object)
    MJDs = COADDS['MJD']
    MJDs[np.where(MJDs.data == '')[0]] = None
 
    COADDS['PROGRAM'] = COADDS['PROGRAM'].astype(object)
    PROGRAM = COADDS['PROGRAM']
    PROGRAM[np.where(PROGRAM.data == '')[0]] = None
    
    COADDS['SDSS_ID'] = COADDS['SDSS_ID'].astype(object)
    CATID = COADDS['SDSS_ID']
    CATID[np.where(CATID.data == '')[0]] = None
 
    if name is not None:
        COADDS = COADDS[np.where(match(COADDS['NAME'].data, name))[0]]
        for coadd in COADDS:
            run1Schema(Table(coadd), topdir, run2d, run1d, DR=coadd['DR'], obs=obs,
                       clobber=clobber, logfile=logfile, use_catid=bool(coadd['USE_CATID']),
                       use_firstcarton=bool(coadd['USE_FIRSTCARTON']),
                       coadd_mjdstart = coadd_mjdstart, usedb=usedb)
    else:
        run1Schema(COADDS, topdir, run2d, run1d, DR=DR, clobber=clobber,
                   logfile=logfile, use_catid=use_catid, obs=obs,
                   use_firstcarton=use_firstcarton, coadd_mjdstart = coadd_mjdstart,
                   usedb=usedb)
    return
    
    
def run1Schema(COADDS, topdir, run2d, run1d, DR=False, clobber=False, logfile=None,
              dronly = False, use_catid=False, use_firstcarton=False, obs=None,
              coadd_mjdstart = None, usedb=False):
    """
    Run all Schema provided by the coaddfile
    """
    if DR:
        for coadd in tqdm(COADDS[np.where(COADDS['DR'] == 1)[0]], desc='DR Coadd',position=0, leave=False, disable=True):
            tqdm.write(str(coadd))
            name = coadd['NAME']
            if obs is not None:
                name = name+'_'+obs
            CustomCoadd(name, topdir, run2d, run1d, cartons=coadd['CARTON'],
                        catalogids=coadd['SDSS_ID'], logfile=logfile, mjd=None,
                        mjdstart=None, mjdend=None, clobber=clobber, obs=obs,
                        rerun1d=coadd['RERUN1D'], use_catid=bool(coadd['USE_CATID']),
                        use_firstcarton=bool(coadd['USE_FIRSTCARTON']),
                        coadd_mjdstart = coadd_mjdstart, usedb=usedb)
    COADDS = COADDS[np.where(COADDS['DR'] == 0)[0]]
    if len(COADDS) == 0:
        splog.info('No Matching Schema')
        return

    topdir2d = ptt.join(topdir, run2d)
    for coadd in tqdm(COADDS, desc='Coadd',position=0, leave=False, disable=True):
        tqdm.write(str(coadd))
        if clean_nones(clean_nones(coadd['MJD'].data)) is None:
            name = coadd['NAME']
            if obs is not None:
                name = name+'_'+obs
            if not clobber:
                frun2ds = glob(ptt.join(field_dir(topdir2d, name, custom=True),
                                        'spPlanCustom-'+name+'-?????.par'))
            else:
                frun2ds = []
            if len(frun2ds) == 0:
                frun2ds = glob(ptt.join(field_dir(topdir2d,'*'), 'spPlan2d-*.par'))
                mjds = [int(ptt.basename(x).split('-')[-1].split('.')[0]) for x in frun2ds]
                if min(mjds) + coadd['CADENCE'] < max(mjds):
                    mjdstart = min(mjds)
                    mjdend = mjdstart+coadd['CADENCE']
            else:
                mjds = [int(ptt.basename(x).split('-')[-1].split('.')[0]) for x in frun2ds]
                mjdstart = max(mjds)+1
                mjdend = mjdstart+coadd['CADENCE']
            if  jdate.astype(int) > mjdend:
                CustomCoadd(name, topdir, run2d, run1d, cartons=coadd['CARTON'], catalogids=coadd['CATID'],
                            logfile=logfile, mjd=None, mjdstart=mjdstart, mjdend=mjdend, clobber=clobber,
                            rerun1d=coadd['RERUN1D'], use_catid=bool(coadd['USE_CATID']), obs=obs,
                            use_firstcarton=bool(coadd['USE_FIRSTCARTON']),
                            coadd_mjdstart = coadd_mjdstart, usedb=usedb)
        else:
            name = coadd['NAME']
            if obs is not None:
                name = name+'_'+obs
            frun2ds = glob(ptt.join(field_dir(topdir2d,name, custom=True),
                                    'spPlanCustom-'+name+'-?????.par'))
            if not clobber:
                mjds = [int(ptt.basename(x).split('-')[-1].split('.')[0]) for x in frun2ds]
            else:
                mjds = []
            if len(mjds) > 0: continue
            
            mjds = clean_nones(clean_nones(coadd['MJD'].data))
            if mjds is None:
                mjds = []
            if len(mjds) > 0:
                if max(mjds) > jdate.astype(int): continue
            CustomCoadd(name, topdir, run2d, run1d, cartons=coadd['CARTON'], catalogids=coadd['CATID'],
                        logfile=logfile, mjd=coadd['MJD'], mjdstart=None, mjdend=None, clobber=clobber,
                        rerun1d=coadd['RERUN1D'], use_catid=bool(coadd['USE_CATID']), obs=obs,
                        use_firstcarton=bool(coadd['USE_FIRSTCARTON']),
                        coadd_mjdstart = coadd_mjdstart, usedb=usedb)
    return
        

########################
