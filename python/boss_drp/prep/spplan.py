#!/usr/bin/env python3
import boss_drp
from boss_drp.field import (field_to_string, Fieldtype, field_dir)
from boss_drp.utils import (find_nearest_indx, Splog, get_dirs, mjd_match, Sphdrfix, getcard)
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
import argparse
try:
    import json
except:
    pass


splog = Splog()

SDSSCOREVersion = getenv('SDSSCORE_VER', default= '')
idlspec2dVersion = boss_drp.__version__

def check_transfer(OBS,mj):
    try:
        evar = f'{OBS.upper()}_STAGING_DATA'
        try:
            transferlog = ptt.join(getenv(evar),'atlogs','{mjd}')
        except:
            transferlog = ptt.join(getenv('APO_STAGING_DATA').replace('apo','{obs}'),
                                            'atlogs','{mjd}')
        transferlog_json = ptt.join(transferlog,'{mjd}_status.json')
        transferlog_json = transferlog_json.format(obs=OBS.lower(), mjd=mj)
        transferlog = ptt.join(transferlog,'transfer-{mjd}.done')
        transferlog = transferlog.format(obs=OBS.lower(), mjd=mj)
        if ptt.exists(transferlog):
            wait = False
        elif ptt.exists(transferlog_json):
            with open(transferlog_json) as lf:
                test = json.load(lf)['history']
                test = {x['stage']:x['status'] for x in test if x['status'] != "skip"}
            if test['copy'] == 'success':
                wait = False
            else:
                wait = True
        else:
            wait = True
    except Exception as e:
        print(e)
        splog.info('Error Checking Data Transfer Logs... continuing anyways')
        wait = False
    return(wait)

def get_alt_cal(fieldexps, allexps, flav='arc', single_cal=False):
    cals  = allexps[np.where(allexps['flavor'].data == flav)[0]].copy()
    idx_n0 = np.where(cals['fieldid'].data != field_to_string(0))[0]

    if len(idx_n0) != 0:
        cals = cals[idx_n0]

    if len(cals) == 0:
        return(fieldexps)
        
    if single_cal is True:
        texp = np.nanmean(fieldexps['TAI'].data)
        idx = find_nearest_indx(cals['TAI'].data, texp)
        cal = cals[idx]
        idx = np.where(cals['EXPOSURE'].data == cal['EXPOSURE'].data)[0]
        cals = cals[idx]
    cals['fieldid'] = fieldexps[0]['fieldid'].data
    fieldexps = vstack([cals,fieldexps])

    return(fieldexps)
    
def get_key(fp):
    filename = ptt.splitext(ptt.splitext(ptt.basename(fp))[0])[0]
    int_part = filename.split('-')[2]
    try:
        return int(int_part)
    except:
        return int(ptt.splitext(int_part)[0])
    

def spplan_findrawdata(inputdir):

    fullnames = glob(ptt.join(inputdir,'sdR*.fit'))
    fullnames.extend(glob(ptt.join(inputdir,'sdR*.fits')))
    fullnames.extend(glob(ptt.join(inputdir,'sdR*.fit.gz')))
    fullnames.extend(glob(ptt.join(inputdir,'sdR*.fits.gz')))
    fullnames = sorted(fullnames, key=get_key)
    fullname_list = []
    keys = []
    for fn in (fullnames):
        key = ptt.basename(fn)
        while '.' in key:
            key = ptt.splitext(key)[0]
        if key in keys:
            continue
        else:
            keys.append(key)
            fullname_list.append(fn)
    return(fullname_list)
        
        





def spplan2d(topdir=None, run2d=None, mjd=None, mjdstart=None, mjdend=None,
             field= None, fieldstart = None, fieldend=None,
             matched_flats=False, nomatched_arcs=False, lco=False, minexp=1,
             clobber=False, release='sdsswork', logfile=None, no_remote=True,
             legacy=False, plates=False, no_commissioning=False, no_dither=False,
             single_flat=False, single_arc=False, override_manual=False, manual_noarc=False,
             verbose = False, returnlist=False, **extra_kwds):
    
    filt_field = field
    if logfile is not None:
        splog = globals()['splog']
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    else:
        if 'splog' in extra_kwds.keys():
            splog = extra_kwds['splog']
        else:
            splog = globals()['splog']
        splog.info('spplan2d started at '+time.ctime())

    if lco:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
        OBS = 'LCO'
        if mjdstart is None:
            mjdstart = 60000
        elif mjdstart <  60000:
            mjdstart = 60000
    else:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_N'
        OBS = 'APO'
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
    
#    sdsscore_dir = getenv('SDSSCORE_DIR')
#    if sdsscore_dir is None:
#        splog.info('ERROR: Must set environment variable SDSSCORE_DIR')
#        exit()
#    sdsscore_dir = ptt.join(sdsscore_dir, OBS.lower())
#    splog.info('Setting SDSSCORE_DIR='+sdsscore_dir)


   #----------
   # Create a list of the MJD directories (as strings)
    mjdlist = get_dirs(rawdata_dir, subdir='', pattern='*', match=mjd, start=mjdstart, end=mjdend)
    nmjd = len(mjdlist)
    splog.info(f'Number of MJDs = {nmjd}')
    
    plateflavors = ['BHM', 'BHM&MWM', 'EBOSS', 'BOSS']
    #---------------------------------------------------------------------------
    # Loop through each input MJD directory

    if returnlist:
        plans_list=[]
    dithered_pmjds = []
    for i, mj in enumerate(mjdlist):
        thismjd = int(mj)

        if OBS == 'APO':
            if thismjd  ==  59560:
                splog.info(f'Skipping {thismjd} for FPS Commissioning')
                continue ##FPS Commissioning
            if thismjd in [59760,59755,59746,59736,59727,59716,59713]:
                splog.info(f'Skipping {thismjd} for 6450Ang Feature')
                continue #6450 Feature:

        ftype = Fieldtype(fieldid=None, mjd=mj)
        if not legacy:
            if ftype.legacy is True:
                continue
        elif not plates:
            if ftype.plates is True:
                continue
        elif not fps:
            if ftype.fps is True:
                continue
        
        inputdir = ptt.join(rawdata_dir, mj)
        sphdrfix = Sphdrfix(mj, fps=ftype.fps, obs=OBS, splog=splog)

#        if ftype.legacy or ftype.plates:
#            plugdir = ptt.join(speclog_dir, mj)
        splog.info('----------------------------')
        splog.info(f'MJD: {mj} {ftype} ({i+1} of {len(mjdlist)})')
        splog.info('Data directory '+inputdir)

        if len(mjdlist) == 1:
            wait = check_transfer(OBS, mj)
            if int(mj) < 59148:
                wait = False
            i = 0
            while wait:
                if i == 0 or i == 1:
                    splog.info('Daily Transfer Log shows incomplete... Waiting 60s')
                else:
                    splog.info('Daily Transfer Log still shows incomplete... continuing anyways')
                    break
                time.sleep(60)
                i = i + 1
                wait = check_transfer(OBS, mj)
                if not wait:
                    break

        # Find all raw FITS files in this directory
        fullname = spplan_findrawdata(inputdir)
        splog.info(f'Number of FITS files found: {len(fullname)}')
        allexps = Table()
        if len(fullname) > 0:
            #-------------------
            # Remove the path from the files names
            
            shortnames = []
            last_confname = ''
            last_conf = None
            for i, f in enumerate(fullname):
                shortnames.append(ptt.basename(f))
            
                #--------------------
                # Find all usefull header keywords
            
                hdr = fits.getheader(f)
                #-----------
                # If mj <= 51813, then set QUALITY=excellent unless this is over-written
                # by SPHDRFIX.  This is because for the early data, this keyword was
                # arbitrarily set to 'unknown', 'bad', 'acceptable', or 'excellent'.
                if thismjd < 51813:
                    hdr.remove('QUALITY', ignore_missing=True, remove_all=True)
                    
                try:
                    expnum = int(ptt.basename(f).split('-')[-1].split('.')[0])
                except:
                    splog.info('ERROR: Cannot determine exposure number from filename '+f)
                    exit(1)
                hdrexp = int(getcard(hdr,'EXPOSURE', default = 0))
                if expnum != hdrexp:
                    splog.info(f'Warning: Exposure number in header({hdrexp}) diagrees with filename ({expnum}) !!')
                    hdr['EXPOSURE'] = expnum

                if not lco:
                    if   (expnum >= 100053) and (expnum <= 100114): hdr['MJD'] = 55050
                    elif (expnum >= 100115) and (expnum <= 100160): hdr['MJD'] = 55051
                    elif (expnum >= 100162) and (expnum <= 100186): hdr['MJD'] = 55052
                    elif (expnum >= 100272) and (expnum <= 100287): hdr['MJD'] = 55061
                    elif (expnum >= 100288) and (expnum <= 100298): hdr['MJD'] = 55062
                    elif (expnum >= 100299) and (expnum <= 100304): hdr['MJD'] = 55063
                    elif (expnum >= 100371) and (expnum <= 100397): hdr['MJD'] = 55065
                    elif (expnum >= 100398) and (expnum <= 100460): hdr['MJD'] = 55066
                    elif (expnum >= 100461) and (expnum <= 100468): hdr['MJD'] = 55067
                    elif (expnum >= 100469) and (expnum <= 100532): hdr['MJD'] = 55068
                    elif (expnum >= 100533) and (expnum <= 100567): hdr['MJD'] = 55069
                    elif (expnum >= 100568) and (expnum <= 100653): hdr['MJD'] = 55070
                    elif (expnum >= 100655) and (expnum <= 100705): hdr['MJD'] = 55071
                    elif (expnum >= 100706) and (expnum <= 100745): hdr['MJD'] = 55072
                    elif (expnum >= 100746) and (expnum <= 100781): hdr['MJD'] = 55073
            
                hdr = sphdrfix.fix(f, hdr)
                FLAVOR = getcard(hdr,'FLAVOR', default='')
                
                #-----------
                # Rename 'target' -> 'science', and 'calibration' -> 'arc'
                try:
                
                    if (FLAVOR.lower() == 'target'):
                        FLAVOR = 'science'
                    elif (FLAVOR.lower() == 'calibration'):
                        FLAVOR = 'arc'
                    elif (FLAVOR.lower() == 'bias'):
                        if verbose:
                            splog.info('Skipping bias file '+ptt.basename(f))
                        continue
                    elif (FLAVOR.lower() == 'dark'):
                        if verbose:
                            splog.info('Skipping dark file '+ptt.basename(f))
                        continue
                    elif (FLAVOR.lower() == ''):
                        if verbose:
                            splog.info('Skipping file '+ptt.basename(f)+' with no flavor')
                        continue
                    if getcard(hdr,'HARTMANN', default='out').lower() in ['right','left']:
                        if (FLAVOR.lower() == 'arc'):
                            if verbose:
                                splog.info('Skipping Hartmann file '+ptt.basename(f))
                            continue
                        else:
                            splog.info(f'Skipping {FLAVOR} file '+ptt.basename(f)+' with closed Hartmann Doors')
                            continue
                except:
                    splog.info(f'Skipping '+ptt.basename(f)+' with failed fits header')
                    continue
                    
                reject = Reject(f, hdr)
                if reject.check(splog):
                    continue
                
                if thismjd > 51576:
                    bin = False
                    if getcard(hdr,'COLBIN', default=1) > 1:
                        splog.info('Skipping binned exposure '+ptt.basename(f))
                        continue
                    elif getcard(hdr,'ROWBIN', default=1) > 1:
                        splog.info('Skipping binned exposure '+ptt.basename(f))
                        continue
                        

                #-----------
                # If the OBSCOMM keyword contains the words "dithered" or "focus",
                # then assume this is test data (set QUALITY=test).
                obscomm = getcard(hdr,'OBSCOMM', default='')
                if ('dithered' in obscomm) or ('focus' in obscomm):
                    hdr['QUALITY'] = 'test'
                
                qual = getcard(hdr,'QUALITY', default='excellent')
                if  (qual != 'excellent') and (qual != "'excellent'"):
                    splog.info('Warning: Non-excellent quality '+FLAVOR+
                              ' file '+ptt.basename(f)+' ('+qual+')')
                    continue
                    
                if ftype.legacy or ftype.plates:
                    dither = 'F'
                    MAPNAME = getcard(hdr,'NAME',default='')
                    fieldid = getcard(hdr,'PLATEID', default=0)
                    platetype = getcard(hdr,'PLATETYP', default='')  #
                    if getcard(hdr,'PLATETYP') is not None:
                        nhdr = hdr.count('PLATETYP')
                    else:
                        nhdr = 0
                    CONFNAME = ''
                    
                    if fieldid == '': fieldid = '0'
                    if str(int(fieldid)) != MAPNAME.split('-')[0]:
                        splog.info('Warning: Plate number '+str(int(fieldid))+
                                  ' flavor '+FLAVOR+ ' inconsistent with map name for' +
                                  ptt.basename(f))
                        continue
                    
                    if len(MAPNAME) <=4:
                        """
                         MAPNAME should be of the form '000000-51683-01'.
                         If it only contains the FIELDID ;; (for MJD <= 51454),
                         then find the actual plug-map file.
                        """
                        plugfile = find_plPlugMapM('*', str(int(MAPNAME)).zfill(4), splog=splog, release=release, no_remote=no_remote)
                        #plugfile = 'plPlugMapM-'+str(int(MAPNAME)).zfill(4)+'-*.par'
                        if plugfile is not None:
                            plugfile = glob(plugfile)
                            #plugfile = glob(ptt.join(plugdir, plugfile))
                            if len(plugfile) == 1:
                                MAPNAME = '-'.join(ptt.basename(plugfile[0]).split('.')[0].split('-')[1:])

                    if (int(fieldid) == 9438) & (thismjd == 58125) & (int(getcard(hdr,'EXPOSURE',default = 0)) == 258988):
                        splog.info('Warning: Skipping Exposure because of trail in data')
                        continue
                else:
                    dither = 'F'
                    MAPNAME = getcard(hdr,'CONFID', default = '0',noNaN=True)
                    fieldid = getcard(hdr,'FIELDID', default = 0, noNaN=True)
                    platetype = 'BHM&MWM'
                    nhdr = 1
                    CONFNAME = MAPNAME

                    if (OBS == 'APO'):
                        if (fieldid in [100520,100542,100981]) and (thismjd == 59733):
                            splog.info(f'Warning: Skipping Exposure because of inadvertent telescope offset')
                            continue
                        if (fieldid == 16165) and (thismjd == 59615):
                            splog.info('Warning: Skipping Exposure because incorrect design/configuration loaded')
                            continue
                        if (fieldid == 20549) and (thismjd == 59623):
                            splog.info('Warning: Skipping Unguided Exposure')
                            continue
                    elif (OBS == 'LCO'):
                        pass

                    if last_confname == CONFNAME:
                        confile = last_conf
                    else:
                        last_confname = CONFNAME
                        confile = get_confSummary(CONFNAME, obs=OBS, splog=splog, release=release,
                                                  sort=False, filter=False, no_remote=no_remote)
                        last_conf = confile
                    if len(confile) == 0:
                        splog.info(f'Warning: Invalid or Missing confSummary for {CONFNAME} for '+str(hdr['EXPOSURE']).strip())
                        if fieldid != 0:
                            continue
                        
                    if 'is_dithered' in confile.meta.keys():
                        if (confile.meta['is_dithered'] == '1') or (confile.meta['parent_configuration'].strip() != '-999'):
                            dither = 'T'
                            if no_dither:
                                splog.info('Warning: Skipping Dither Exposure '+str(hdr['EXPOSURE']).strip())
                                continue
                        
                    if 'field_id' in confile.meta.keys():
                        if (confile.meta['field_id'].strip() == '-999'):
                            fieldid = field_to_string(0)
                        if len(confile.meta['field_id'].strip()) == 0:
                            fieldid = field_to_string(0)
                    
                    if FLAVOR == 'science':
                        ftype_exp = Fieldtype(fieldid=fieldid, mjd=mj)
                        if (ftype_exp.engineering):
                            splog.info('Warning: Skipping Engineering Exposure '+str(hdr['EXPOSURE']).strip())
                            continue
                    
                if platetype not in plateflavors:
                    splog.info('Skipping '+platetype+' field '+str(fieldid)+' exposure '+str(hdr['EXPOSURE']).strip())
                    continue
                    
                if FLAVOR == 'science':
                    try:
                        ffs = np.asarray(getcard(hdr,'FFS', default = '1 1 1 1 1 1 1 1').split(), dtype=int)
                    except:
                        splog.info(f'Warning: {ptt.basename(f)} with undetermined Flat Field Sutter Status')
                        ffs = np.asarray([0,0,0,0,0,0,0,0])
                    if sum(ffs) != 0:
                        splog.info('Warning: Flat Field Shutters closed for science exposure '+ptt.basename(f))
                        continue
                elif (FLAVOR == 'arc') or (FLAVOR == 'flat'):
                    try:
                        ffs = np.asarray(getcard(hdr,'FFS', default = '0 0 0 0 0 0 0 0').split(), dtype=int)
                    except:
                        splog.info(f'Warning: {ptt.basename(f)} with undetermined Flat Field Sutter Status')
                        ffs = np.asarray([1,1,1,1,1,1,1,1])
                    if sum(ffs) != len(ffs):
                        splog.info(f'Warning: Flat Field Shutters open for {FLAVOR} exposure '+ptt.basename(f))
                        continue

                if int(getcard(hdr,'MJD', default = 0)) != thismjd:
                    splog.info('Warning: Wrong MJD in file '+ptt.basename(f))
                
                if FLAVOR.upper() == 'UNKNOWN':
                    continue
                
                if ftype.legacy or ftype.plates:
                    fid = str(fieldid).strip()
                else:
                    fid = field_to_string(fieldid)
                exp = Table({'exptime':[np.float32(getcard(hdr,'EXPTIME', default = 0.0))],
                            'EXPOSURE':[int(getcard(hdr,'EXPOSURE', default=0))],
                            'flavor':[FLAVOR],
                            'DITHER':[dither],
                            'mapname':[str(MAPNAME)],
                            'TAI':[float(getcard(hdr, 'TAI-BEG', default = 0.0))],
                            'mjd':[np.int32(thismjd)],
                            'fieldid':[fid],
                            'confid':[str(CONFNAME)],
                            'CONFNAME':[str(CONFNAME)],
                            'shortname':[ptt.basename(f)]})

                allexps = vstack([allexps,exp])
        if len(allexps) == 0:
            splog.info('No Science Frames for '+mj)
        else:
            if ftype.legacy or ftype.plates:
                fieldmap_col = 'mapname'
            else:
                fieldmap_col = 'fieldid'
            for field in list(dict.fromkeys(allexps[fieldmap_col].data)):
            
                manual = 'F'
                if ftype.fps:
                    if int(field) == 0:
                        continue
                if ftype.legacy or ftype.plates:
                    ftest = field.split('-')[0]
                else:
                    ftest = field
                if not mjd_match(ftest, mjd=filt_field, mjdstart=fieldstart, mjdend=fieldend):
                    splog.info(f'Skipping Field {ftest} outside specified field range/list')
                    continue
                # Filter to a single FPS Field or Plate Map
                fieldexps = allexps[np.where(allexps[fieldmap_col].data == field)[0]]
                sci = (np.logical_or((np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'science'),
                                     (np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'smear')))
                nsci = len(fieldexps[np.where(sci)[0]])
                if nsci == 0:
                    # Check for valid science frames
                    splog.info(f'WARNING: No science frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue
                elif nsci < minexp:
                    splog.info(f'WARNING: Insufficient ({nsci}<{minexp}) science frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue

                if not (ftype.legacy or ftype.plates):
                    if len(fieldexps[np.where((fieldexps['flavor'].data == 'arc'))[0]]) == 0:
                        # Check for valid arc Frame
                        if nomatched_arcs:
                            fieldexps = get_alt_cal(fieldexps, allexps, flav='arc', single_cal=single_arc)
                        elif manual_noarc:
                            manual = 'T'
                            pf = f'spPlan2d-{field_to_string(field)}-{mj}.par'
                            splog.info(f'WARNING: Building plan {pf} for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd}) as manual with unmatched arcs')
                            fieldexps = get_alt_cal(fieldexps, allexps, flav='arc', single_cal=single_arc)
                if len(fieldexps[np.where((fieldexps['flavor'].data == 'arc'))[0]]) == 0:
                    splog.info(f'WARNING: No arc frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue


                if not (ftype.legacy or ftype.plates):
                    if len(fieldexps[np.where((fieldexps['flavor'].data == 'flat'))[0]]) == 0:
                        # Check for valid flat Frame
                        if not matched_flats:
                            fieldexps = get_alt_cal(fieldexps, allexps, flav='flat', single_cal=single_arc)
                if len(fieldexps[np.where((fieldexps['flavor'].data == 'flat'))[0]]) == 0:
                    splog.info(f'WARNING: No flat frames for {fieldmap_col} {field_to_string(field)} (mjd:{thismjd})')
                    continue

                sci = (np.logical_or((np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'science'),
                                     (np.char.strip(np.asarray(fieldexps['flavor'].data)) == 'smear')))
                fieldname = field_to_string(fieldexps[sci]['fieldid'].data[0])
                
                ftype_exp = Fieldtype(fieldid=fieldname, mjd=mj)
                if ftype.legacy is not ftype_exp.legacy:
                    splog.info('Warning: Skipping Legacy plate from non-Legacy MJD')
                    continue
                elif ftype.plates is not ftype_exp.plates:
                    splog.info('Warning: Skipping SDSS-V plate from non-SDSS-V Plate MJD')
                    continue
                elif ftype_exp.commissioning is True and no_commissioning is True:
                    splog.info('Warning: Skipping SDSS-V FPS commissioning Field')
                    continue
                elif ftype.fps is not ftype_exp.fps:
                    splog.info('Warning: Skipping FPS Field from non-FPS MJD')
                    continue
                DITHER = fieldexps[sci]['DITHER'].data[0]
                planfile = 'spPlan2d-' + fieldname + '-' + mj + '.par'
                planfile = ptt.join(field_dir(ptt.join(topdir,run2d), fieldname), planfile)

                if returnlist:
                    plans_list.append(planfile)
                if ftype.legacy:
                    shape = (4,)
                else:
                    shape = (2,)
                
                fieldexps.add_column(Column('', dtype=object, name='name', shape=shape))

                for i, row in enumerate(fieldexps):
                    tnames = fieldexps[np.where(fieldexps['EXPOSURE'].data == row['EXPOSURE'])[0]]['shortname'].data
                    tnames = np.char.replace(tnames, '.gz', '')
                    names = []
                    for cam in ['b1', 'b2', 'r1', 'r2']:
                        for tn in tnames:
                            if cam in tn:
                                names.append(tn)
                    while len(names) < shape[0]:
                        names.append('')
                    fieldexps[i]['name'] = np.asarray(names)
                names = fieldexps['name'].data
                names = np.stack(names.tolist(), axis=0)
                names = names.astype(str)
                fieldexps.remove_column('name')
                fieldexps.add_column(Column(names, dtype=str, name='name', shape=shape))
                fieldexps = unique(fieldexps, keys='EXPOSURE')

                fieldexps = fieldexps['confid','fieldid','mjd','mapname','flavor','exptime','name']
        
                
                fieldexps.meta=OrderedDict({
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
                
                if ptt.exists(planfile):
                    if clobber is False:
                        splog.info('WARNING: Will not over-write plan file: ' + ptt.basename(planfile))
                        continue
                    else:
                        test = read_table_yanny(planfile, 'SPEXP')
                        if not override_manual:
                            if 'manual' in test.meta.keys():
                                if test.meta['manual'] == 'T':
                                    splog.info('WARNING: Will not over-write manual plan file: ' + ptt.basename(planfile))
                                    continue
                        else:
                            if 'manual' in test.meta.keys():
                                if test.meta['manual'] == 'T':
                                    splog.info('Backing up manual plan file to '+ptt.splitext(planfile)[0]+'.bkup')
                                    rename(planfile, ptt.splitext(planfile)[0]+'.bkup')
                        splog.info('WARNING: Over-writing plan file: ' + ptt.basename(planfile))
                else:
                    splog.info('Writing plan file '+ ptt.basename(planfile))
                makedirs(ptt.dirname(planfile), exist_ok = True)
                fieldexps.convert_unicode_to_bytestring()
                write_table_yanny(fieldexps, planfile,tablename='SPEXP', overwrite=clobber)
                del fieldexps
        del allexps
    splog.info('----------------------------')
    splog.info('Successful completion of spplan2d at '+ time.ctime())
    if 'skip1d' not in extra_kwds.keys():
        extra_kwds['skip1d'] = False

    if extra_kwds['skip1d'] is True and logfile is not None:
        splog.close()
    if returnlist:
        return(plans_list)
    else:
        return

def spplan1d (topdir=None, run2d=None, mjd=None, mjdstart=None, mjdend=None,
             field= None, fieldstart = None, fieldend=None, lco=False,
             clobber=False, logfile=None, override_manual=False,
             legacy=False, plates=False, plate_epoch = False,
             daily = False, plans=None, **extra_kwds):
    
    if plate_epoch is False: daily=True
    if logfile is not None and extra_kwds['skip2d'] is True:
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    else:
        if 'splog' in extra_kwds.keys():
            splog = extra_kwds['splog']
        else:
            splog = globals()['splog']
        splog.info('spplan1d started at '+time.ctime())
    
    #----------
    # Determine the top-level of the directory tree
    #----------
    if topdir is None:
        topdir = getenv('BOSS_SPECTRO_REDUX')
    splog.info('Setting TOPDIR='+ topdir)

    if run2d is None:
           run2d = getenv('RUN2D')
    splog.info('Setting RUN2D='+ run2d)
    topdir2d = ptt.join(topdir, run2d)

    if not(ptt.exists(topdir) and ptt.isdir(topdir)):
        splog.info('Directory does not exist: '+topdir)
        exit()
   
    
    OBS = 'LCO' if lco else 'APO'
    
    if plans is not None:
        if field is None:
            field = []
        field.extend([ptt.basename(x).split('-')[1] for x in np.atleast_1d(plans)])
    

    fieldlist = get_dirs(ptt.dirname(field_dir(topdir2d, '*')), field = True,
                         match=field, start=fieldstart, end=fieldend)
    splog.info('Number of field directories = '+ str(len(fieldlist)))

    # Loop through each input configuration directory
    for fielddir in fieldlist:
        try: 
            fieldid = int(ptt.basename(fielddir))
        except:
            continue
        ftype = Fieldtype(fieldid=fieldid, mjd=mjd)
        splog.info('----------------------------')
        splog.info('Field directory '+field_dir(topdir2d, fielddir))
        #----------
        # Find all 2D plan files
        allplan = glob(ptt.join(field_dir(topdir2d, fielddir), 'spPlan2d*.par'))
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
                test = mjd_match(row['mjd'], mjd=mjd, mjdstart=mjdstart, mjdend=mjdend)
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
                    planfile = ptt.join(field_dir(topdir2d, fielddir), planfile)

                    fmjds_exps.meta=OrderedDict({
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
                                'SDSS_access_Ver':  "'"+saver+"'"             +"   # Version of sdss_access when building plan file",
                                'manual':           "F"                       +"   # Manually edited plan file (T: True, F: False)"
                                         })
                
                    if ptt.exists(planfile):
                        if clobber is False:
                            splog.info('WARNING: Will not over-write plan file: ' + ptt.basename(planfile))
                            continue
                        else:
                            test = read_table_yanny(planfile, 'SPEXP')
                            if not override_manual:
                                if 'manual' in test.meta.keys():
                                    if test.meta['manual'] == 'T':
                                        splog.info('WARNING: Will not over-write manual plan file: ' + ptt.basename(planfile))
                                        continue
                            else:
                                if 'manual' in test.meta.keys():
                                    if test.meta['manual'] == 'T':
                                        splog.info('Backing up manual plan file to '+ptt.splitext(planfile)[0]+'.bkup')
                                        rename(planfile, ptt.splitext(planfile)[0]+'.bkup')
                            splog.info('WARNING: Over-writing plan file: ' + ptt.basename(planfile))
                    else:
                        splog.info('Writing plan file '+ ptt.basename(planfile))
                    fmjds_exps.convert_unicode_to_bytestring()
                    write_table_yanny(fmjds_exps, planfile,tablename='SPEXP', overwrite=clobber)
    splog.info('----------------------------')
    splog.info('Successful completion of spplan1d at '+ time.ctime())
    if logfile is not None:
        splog.close()
    return
                
                
