#!/usr/bin/env python3
from boss_drp import __version__ as idlspec2dVersion
from boss_drp.utils.splog import splog
from boss_drp.prep.spplan.tools import *
from boss_drp.field import (field_to_string, Fieldtype, Field)
from boss_drp.utils import (get_dirs, mjd_match)
from boss_drp.prep import check_manual_cal
from boss_drp.Config import config

from sdss_access import __version__ as saver
from tree import __version__ as treever

from os import getenv
import os.path as ptt
from glob import glob
import time
from astropy.table import Table, vstack, Column
from collections import OrderedDict
from pydl.pydlutils.yanny import read_table_yanny
from pydl import __version__ as pydlVersion
import numpy as np

SDSSCOREVersion = getenv('SDSSCORE_VER', default= '')

def spplan2d(pipe=False):
    
    logfile = config.pipe['plan.daily.dailyplan_logfile']
    if logfile is not None:
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    splog.info('spplan2d started at '+time.ctime())

    filt_field = config.pipe['fmjdselect.field']
    fieldrange = config.pipe['fmjdselect.fieldrange']

    if type(config.pipe['fmjdselect.obs']) == list:
        config.pipe['fmjdselect.obs'] = config.pipe['fmjdselect.obs'][0]

    mjdrange = config.pipe['fmjdselect.mjdrange']
    if mjdrange is None:
        mjdrange = [[None, None]]
    if not isinstance(mjdrange[0], list):
        mjdrange = [mjdrange]
    if config.pipe['fmjdselect.obs'].lower() == 'lco':
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
        OBS = 'LCO'
        for i, mjdr in enumerate(mjdrange):
            if mjdr[0] is None:
                mjdrange[i][0] = 60000
            elif int(mjdr[0]) < 60000:
                mjdrange[i][0] = 60000
            else:
                mjdrange[i][0] = int(mjdr[0])
            if mjdr[1] is not None:
                mjdrange[i][1] = int(mjdr[1])
        config.pipe['fmjdselect.mjdrange'] = mjdrange
    else:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_N'
        OBS = 'APO'


    lco = False if OBS == 'APO' else True
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
    mjdlist = get_dirs(rawdata_dir, subdir='', pattern='*', 
                       match=config.pipe['fmjdselect.mjd'], 
                       ranges = config.pipe['fmjdselect.mjdrange'])
    nmjd = len(mjdlist)
    splog.info(f'Number of MJDs = {nmjd}')
    if nmjd == 0:
        splog.info('No Valid MJDs')
        return None
    #---------------------------------------------------------------------------
    # Loop through each input MJD directory

    plans_list=[]
    dithered_pmjds = []

    legacy = config.pipe['SDSS_Generation.legacy']
    plates = config.pipe['SDSS_Generation.plates']
    fps = config.pipe['SDSS_Generation.fps']
    if config.pipe['SDSS_Generation.sdssv'] is True:
        fps = True
        plates = True
    no_dither = not config.pipe['fmjdselect.dither']
    no_remote = not config.pipe['general.REMOTE']
    release = config.pipe['general.RELEASE']
    verbose = config.pipe['plan.daily.daily_plan_verbose']
    manual_noarc = config.pipe['plan.daily.flag_nomatch_manual']
    minexp = config.pipe['plan.daily.minscience']
    nomatched_arcs =  not config.pipe['plan.daily.matched_arc']
    matched_flats = config.pipe['plan.daily.matched_flat']
    single_flat = not config.pipe['plan.daily.multiple_flat']
    single_arc = not config.pipe['plan.daily.multiple_arcs']

    clobber_traceplan = config.pipe['Clobber.clobber_spTrace']
    run_write_traceplan = config.pipe['plan.daily.run_traceplan'] if not pipe else False
    
    override_manual_trace = config.pipe['pipe.trace.override_manual_trace']
    for i, mj in enumerate(mjdlist):
        ftype = Fieldtype(fieldid=None, mjd=mj, obs=OBS)
        if not ftype.check(legacy=legacy, plates=plates, fps=fps):
            splog.info(f'Skipping {mj}: {ftype}')
            continue
        splog.info('----------------------------')
        splog.info(f'MJD: {mj} {ftype} ({i+1} of {len(mjdlist)})')

        allexps, ftype = build_exps(i, mj, mjdlist, OBS, rawdata_dir, ftype, spplan_Trace=False,
                                    legacy=legacy, plates=plates, fps=fps,
                                    lco=lco,no_dither=no_dither,
                                    no_remote=no_remote, release=release, verbose=verbose)
        
        thismjd = int(mj)
        if len(allexps) == 0:
            splog.info('No Science Frames for '+mj)
        else:
            if ftype.legacy or ftype.plates:
                fieldmap_col = 'mapname'
            else:
                fieldmap_col = 'fieldid'
            manual_noarc_set = manual_noarc
            
            traceflat = False
            if lco:
                if int(thismjd) >= obsTrace_mjdstart['LCO']:
                    traceflat = True
            else:
                if int(thismjd) >= obsTrace_mjdstart['APO']:
                    traceflat = True


            if run_write_traceplan and traceflat:
                allexps_trace = allexps[np.isin(allexps["flavor"], ['arc','flat','calibration'])]

                write_traceplan(allexps_trace, mj, ftype, OBS, False, None, topdir, 
                        run2d, release, clobber_traceplan, override_manual_trace, plan2d=True)


                del allexps_trace


            if traceflat:
                allexps = get_master_cal(allexps, dropMaster = False, obs= OBS, mjd=mj)
                if allexps is None:
                    continue
            for field in list(dict.fromkeys(allexps[fieldmap_col].data)):
                if ftype.legacy or ftype.plates:
                    ftype_exp = Fieldtype(mjd=mj, obs='apo')
                else:
                    ftype_exp = Fieldtype(fieldid=field_to_string(field), mjd=mj, obs=OBS)
                if ftype.legacy is not ftype_exp.legacy:
                    splog.info(f'Warning: Skipping Legacy plate {field_to_string(field)} from non-Legacy MJD')
                    continue
                elif ftype.plates is not ftype_exp.plates:
                    splog.info(f'Warning: Skipping SDSS-V plate {field_to_string(field)} from non-SDSS-V Plate MJD')
                    continue
                elif ftype_exp.commissioning is True and config.pipe['fmjdselect.commissioning'] is False:
                    splog.info(f'Warning: Skipping SDSS-V FPS commissioning Field {field_to_string(field)}')
                    continue
                elif ftype.fps is not ftype_exp.fps:
                    splog.info(f'Warning: Skipping FPS Field {field_to_string(field)} from non-FPS MJD')
                    continue

                manual_noarc = manual_noarc_set
                ## The code below handles a small number of cases where a field does not have a valid arc,
                ## and uses the same as another field from the same night

                manual_noarc, use_arc = check_manual_cal(type='arc', field=field,
                                                         mjd=thismjd,
                                                         obs = 'lco' if lco else 'apo')
                manual_noflag, use_flat = check_manual_cal(type='flat', field=field,
                                                         mjd=thismjd,
                                                         obs = 'lco' if lco else 'apo')
                manual = 'F'
                if ftype.fps:
                    if int(field) == 0:
                        continue
                if ftype.legacy or ftype.plates:
                    ftest = field.split('-')[0]
                else:
                    ftest = field
                if not mjd_match(ftest, mjd=filt_field, ranges = fieldrange):
                    splog.info(f'Skipping Field {ftest} outside specified field range/list')
                    continue
                # Filter to a single FPS Field or Plate Map
                fieldexps = allexps[np.where(allexps[fieldmap_col].data == field)[0]]
                if use_arc is not None:
                    mask = (fieldexps['flavor'] != 'arc') | (fieldexps['EXPOSURE'] == use_arc)
                    fieldexps = fieldexps[mask]
                if use_flat is not None:
                    mask = (fieldexps['flavor'] != 'flat') | (fieldexps['EXPOSURE'] == use_flat)
                    fieldexps = fieldexps[mask]
                    
                sci = (np.logical_or((np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'science'),
                                     (np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'smear')))
                nsci = len(fieldexps[np.where(sci)[0]])
                if nsci == 0:
                    # Check for valid science frames
                    splog.info(f'WARNING: No science frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue
                elif nsci < (minexp or 1):
                    splog.info(f'WARNING: Insufficient ({nsci}<{minexp}) science frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue
                splog.info(f'Building Plan for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                fieldexps = check_cal_dt(fieldexps)
                if lco:
                    fieldexps = check_cal_Screen(fieldexps, flavor='arc')
                    fieldexps = check_cal_Screen(fieldexps, flavor='flat')
                if traceflat:
                    fieldexps = get_alt_cal(fieldexps, allexps, flav='TRACEFLAT')
                    fieldexps = get_alt_cal(fieldexps, allexps, flav='TRACEARC')

                if not (ftype.legacy or ftype.plates):
                    if len(fieldexps[np.where((fieldexps['flavor'].data == 'arc'))[0]]) == 0:
                        # Check for valid arc Frame
                        if nomatched_arcs:
                            fieldexps = get_alt_cal(fieldexps, allexps, flav='arc', single_cal=single_arc)
                        elif manual_noarc:
                            manual = 'T'
                            pf = f'spPlan2d-{field_to_string(field)}-{mj}.par'
                            msg = f'WARNING: Building plan {pf} for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd}) as manual with unmatched arcs'
                            fieldexps = get_alt_cal(fieldexps, allexps, flav='arc', single_cal=single_arc, use_cal=use_arc, msg= msg)
                if len(fieldexps[np.where((fieldexps['flavor'].data == 'arc'))[0]]) == 0:
                    splog.info(f'WARNING: No arc frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue


                if not (ftype.legacy or ftype.plates):
                    if len(fieldexps[np.where((fieldexps['flavor'].data == 'flat'))[0]]) == 0:
                        # Check for valid flat Frame
                        if not matched_flats:
                            fieldexps = get_alt_cal(fieldexps, allexps, flav='flat', single_cal=single_flat, use_cal = use_flat)
                if len(fieldexps[np.where((fieldexps['flavor'].data == 'flat'))[0]]) == 0:
                    splog.info(f'WARNING: No flat frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue

                sci = (np.logical_or((np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'science'),
                                     (np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'smear')))
                fieldname = field_to_string(fieldexps[sci]['fieldid'].data[0])
                
                DITHER = fieldexps[sci]['DITHER'].data[0]
                planfile = 'spPlan2d-' + fieldname + '-' + mj + '.par'
                fc = Field(topdir, run2d, fieldname)
                planfile = ptt.join(fc.dir(), planfile)

                plans_list.append(planfile)
                meta = OrderedDict({
                            'fieldname':        fieldname                +"   # Field number",
                            'MJD':              mj                       +"   # Modified Julian Date",
                            'OBS':              OBS                      +"   # Observatory",
                            'RUN2D':            run2d                    +"   # 2D reduction name",
                            'DITHER':           DITHER                   +"   # Is the Field Dithered (T: True, F: False)",
                            'planfile2d': "'"+ptt.basename(planfile)+"'" +"   # Plan file for 2D spectral reductions (this file)",
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
                fieldexps = pair_ccds(ftype, fieldexps, OBS=OBS)
                write_plan(planfile, fieldexps, meta=meta, clobber=config.pipe['Clobber.clobber_plan'], 
                           override_manual=config.pipe['plan.daily.override_manual'])
                del fieldexps
        del allexps
    splog.info('----------------------------')
    splog.info('Successful completion of spplan2d at '+ time.ctime())

    if config.pipe['plan.daily.skip1d']:
        if logfile is not None:
            splog.close()
    elif config.pipe['plan.daily.quick1d']:
        return(plans_list)
    return

def spplan1d (plans):
    
    if not config.pipe['plan.daily.plate_epoch']:
        daily = True

    logfile = config.pipe['plan.daily.dailyplan_logfile']
    if (config.pipe['plan.daily.skip2d']) and (logfile is not None):
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    splog.info('spplan1d started at '+time.ctime())
    
    #----------
    # Determine the top-level of the directory tree
    #----------
    topdir = config.pipe['general.BOSS_SPECTRO_REDUX']
    splog.info('Setting TOPDIR='+ topdir)

    run2d = config.pipe['general.RUN2D']
    splog.info('Setting RUN2D='+ run2d)

    if not(ptt.exists(topdir) and ptt.isdir(topdir)):
        splog.info('Directory does not exist: '+topdir)
        exit()
   
    
    OBS = config.pipe['fmjdselect.obs'].upper()

    field = None
    frange = None
    if plans is not None:
        field = config.pipe['fmjdselect.field']
        if field is None:
            field = []
        field.extend([ptt.basename(x).split('-')[1] for x in np.atleast_1d(plans)])
        frange = config.pipe['fmjdselect.fieldrange']

    mjd = np.atleast_1d(config.pipe['fmjdselect.mjd'])
    afc = Field(topdir, run2d, '*')
    fieldlist = get_dirs(ptt.dirname(afc.dir()), field = True,
                         match=field, ranges = frange)
    splog.info('Number of field directories = '+ str(len(fieldlist)))

    # Loop through each input configuration directory
    for fielddir in fieldlist:
        try: 
            fieldid = int(ptt.basename(fielddir))
        except:
            continue
        fc = Field(topdir, run2d, fielddir, mjd=mjd, obs=OBS)
        ftype = fc.type
        splog.info('----------------------------')
        splog.info('Field directory '+fc.dir())
        #----------
        # Find all 2D plan files
        allplan = glob(ptt.join(fc.dir(), 'spPlan2d*.par'))
        #----------
        # Read all the 2D plan files
        # The string array PLANLIST keeps a list of the plan file that each element
        # of the ALLEXP structure came from, and MJDLIST keeps the list of each MJD
        allexp = Table()
        for thisplan in allplan:
            thisexp = read_table_yanny(thisplan, 'SPEXP')
            thisexp.convert_bytestring_to_unicode()
            
            try:
                if thisexp.meta['OBS'] != OBS:
                    continue
            except:
                if thisexp['name'][0][0].split('-')[1] in ['b2','r2']:
                    tobs = 'LCO'
                else:
                    tobs = 'APO'
                if tobs != OBS:
                    continue

            sci = (np.logical_or((np.char.strip(np.asarray(thisexp['flavor'].data)) == 'science'),
                                 (np.char.strip(np.asarray(thisexp['flavor'].data)) == 'smear')))
            thisexp = thisexp[np.where(sci)[0]]
            try:
                thisexp.add_column(thisexp.meta['DITHER'], name='DITHER')
            except:
                thisexp.add_column('F', name='DITHER')
            thisexp.add_column(Column(ptt.basename(thisplan),name='thisplan', dtype=object))
            thisexp.meta = {}
            allexp = vstack([allexp, thisexp])
        if len(allexp) == 0:
            splog.info(f'No valid plans for {OBS}')
            continue
        if ftype.legacy or ftype.plates:
            fieldmap_col = 'mapname'
        else:
            fieldmap_col = 'fieldid'
        for fld in list(dict.fromkeys(allexp[fieldmap_col].data)):
            # Filter to a single FPS Field or Plate Map
            spexp = allexp[np.where(allexp[fieldmap_col].data == fld)[0]]
            if len(spexp) == 0:
                continue
            badMjds=[]
            # Decide if any of these MJD's are within the bounds specified by MJD,MJSTART,MJEND.
            for i,row in enumerate(spexp):
                test = mjd_match(row['mjd'], mjd=config.pipe['fmjdselect.mjd'], ranges = config.pipe['fmjdselect.mjdrange'])
                if test is False:
                    badMjds.append(i)
            if len(badMjds) > 0:
                spexp.remove_rows(badMjds)
            if len(spexp) > 0:
                # -------------
                # Replace the prefix 'sdR' with 'spFrame' in the science frames
                # and the suffix '.fit' with '.fits'
                names = spexp['name'].data
                names = np.char.replace(names, 'sdR', 'spFrame')
                names = np.char.replace(names, '.fit', '.fits')
                spexp.remove_column('name')
                spexp.add_column(names, name='name')
                
                if daily:
                    epoch_len = 1
                else:
                    if ftype.rm_plate:
                        epoch_len = 3
                    elif ftype.plates or ftype.legacy:
                        epoch_len = 1000
                    else:
                        epoch_len = 1
                spexp.add_column(np.int32(-1), name='epoch_combine')
                
                while len(np.where((spexp['epoch_combine'].data == -1))[0]) > 0:
                    nomatch_idx = np.where((spexp['epoch_combine'].data == -1))[0]
                    epoch = np.min(spexp[nomatch_idx]['mjd'].data)
                    idx = np.where((spexp['mjd'].data < epoch+epoch_len) &
                                   (spexp['epoch_combine'].data == -1))[0]
                    ec = spexp['epoch_combine']
                    ec[idx] = np.int32(epoch)
                    try:
                        DITHER = 'T' if 'T' in (spexp[idx]['DITHER'].data) else 'F'
                    except:
                        DITHER = 'F'
                    plan2dfiles = "'"+"' '".join(np.unique(spexp[idx]['thisplan'].data).tolist())+"'"
                    fmjds_exps = spexp[idx]['confid','fieldid','mjd','mapname','flavor','exptime','name', 'epoch_combine']
                    coadd_mjd = np.max(fmjds_exps['mjd'].data)
                    planfile = 'spPlancomb-' + field_to_string(fieldid) + '-' + str(coadd_mjd) + '.par'
                    planfile = ptt.join(fc.dir(), planfile)

                    meta = OrderedDict({
                                'fieldid': field_to_string(fieldid)           +"   # Field number",
                                'MJD':              str(coadd_mjd)            +"   # Modified Julian Date",
                                'OBS':              OBS                       +"   # Observatory",
                                'RUN2D':            run2d                     +"   # 2D reduction name",
                                'DITHER':           DITHER                    +"   # Is the Field Dithered (T: True, F: False)",
                                'planfile2d':       plan2dfiles               +"   # Plan file for 2D spectral reductions",
                                'planfilecomb':"'"+ptt.basename(planfile)+"'" +"   # Plan file for coadding (this file)",
                                'idlspec2dVersion': "'"+idlspec2dVersion+"'"  +"   # Version of idlspec2d when building plan file",
                                #'idlutilsVersion':  "'"+idlutilsVersion+"'"   +"   # Version of idlutils when building plan file",
                                'pydlVersion':      "'"+pydlVersion+"'"       +"   # Version of pydl when building plan file",
                                #'speclogVersion':   "'"+speclogVersion+"'"    +"   # Version of speclog when building plan file",
                                'SDSSCOREVersion':  "'"+SDSSCOREVersion+"'"   +"   # Version of SDSSCORE when building plan file",
                                'SDSS_access_Ver':  "'"+saver+"'"             +"   # Version of sdss_access when building plan file",
                                'manual':           "F"                       +"   # Manually edited plan file (T: True, F: False)"
                                         })
                    write_plan(planfile, fmjds_exps, meta=meta, clobber=config.pipe['Clobber.clobber_plan'], 
                               override_manual=config.pipe['plan.daily.override_manual'])

    splog.info('----------------------------')
    splog.info('Successful completion of spplan1d at '+ time.ctime())
    if logfile is not None:
        splog.close()
    return
                
                
