#!/usr/bin/env python3
import boss_drp
from boss_drp.utils.splog import splog
from boss_drp.field import (field_to_string, Fieldtype, Field)
from boss_drp.utils import (find_nearest_indx, get_dirs, mjd_match, Sphdrfix, getcard)
from boss_drp.prep.GetconfSummary import find_confSummary, find_plPlugMapM, get_confSummary
from boss_drp.utils.reject import Reject
from boss_drp.prep import check_manual_cal
from boss_drp.Config import config

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



SDSSCOREVersion = getenv('SDSSCORE_VER', default= '')
idlspec2dVersion = boss_drp.__version__
obsTrace_mjdstart = {'LCO':60187, 'APO':59560}
plateflavors = ['BHM', 'BHM&MWM', 'EBOSS', 'BOSS']

def check_transfer(OBS,mj):
    """ Check if the Observatory to Utah Transfer is complete """
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

def get_alt_cal(fieldexps, allexps, flav='arc', single_cal=False, use_cal = None, msg = None):
    """ Get alternative Calibration frames that are no natively associated with the field"""
    if msg is not None:
        splog.info(msg)
    cals  = allexps[np.where(allexps['flavor'].data == flav)[0]].copy()
    if use_cal is None:
        idx_n0 = np.where(cals['fieldid'].data != field_to_string(0))[0]

        if len(idx_n0) != 0:
            cals = cals[idx_n0]

    if len(cals) == 0:
        return(fieldexps)
    
    if use_cal is not None:
        splog.info(f'Using ExposureID={use_cal} as {flav.upper()}')
        cals = cals[cals['EXPOSURE'] == use_cal]
        #cals['fieldid'] = fieldexps[0]['fieldid'].data
        fieldexps = vstack([cals,fieldexps])
        return(fieldexps)
        
    if single_cal is True:
        texp = np.nanmean(fieldexps['TAI'].data)
        idx = find_nearest_indx(cals['TAI'].data, texp)
        cal = cals[idx]
        idx = np.where(cals['EXPOSURE'].data == cal['EXPOSURE'].data)[0]
        cals = cals[idx]
    #cals['fieldid'] = fieldexps[0]['fieldid'].data
    for c in cals:
        if str(c['confid']) in ['0','-999']:
            c['confid'] = fieldexps[0]['confid']
        if str(c['mapname']) in ['0','-999']:
            c['mapname'] = fieldexps[0]['mapname']
    fieldexps = vstack([cals,fieldexps])

    return(fieldexps)
    
def get_key(fp):
    """ Returns the exposure ID portion of the filename for sorting key"""
    filename = ptt.splitext(ptt.splitext(ptt.basename(fp))[0])[0]
    int_part = filename.split('-')[2]
    try:
        return int(int_part)
    except:
        return int(ptt.splitext(int_part)[0])
    
    
def check_cal_dt(fieldexps, dt = (2*60*60)):#2hrs)
    """ Check if Flats are more then dt seconds away from the arc and science frames """
    flat     = fieldexps[fieldexps['flavor'] == 'flat']
    if len(flat) == 0:
        return fieldexps
    arc      = fieldexps[fieldexps['flavor'] == 'arc']
    notFlat  = fieldexps[fieldexps['flavor'] != 'flat']

    field = arc if len(arc) > 0 else notFlat
    differences = np.abs(flat['TAI'].data[:, None] - field['TAI'].data)

    # Check if any difference is less than dt
    if np.any(differences <= dt):
        splog.debug(f'Minumum dTau = {round(np.min(differences),2)}s')
        return fieldexps
    splog.info(f'Dropping associated flats... All flats > {dt}s from Field (Sci/Arc) exposures')
    splog.debug(f'Minumum dTau = {round(np.min(differences),2)}s')
    return notFlat


def mark_LCOfixedScreenCals(fieldexps, limit = 40.0):
    """ Mark Frames taken with the LCO Fixed Dome Flat Screen
        The windscreen has a lower Alt lim, so any cals below
        the limit must be with the fixed screen """
    if 'FixedScreen' in fieldexps.columns:
        return fieldexps
    fieldexps['FixedScreen'] = False
    mask = fieldexps['alt'] < float(limit)
    fieldexps['FixedScreen'][mask] = True
    return fieldexps
    
def check_cal_Screen(fieldexps, flavor= 'arc'):
    fieldexps = mark_LCOfixedScreenCals(fieldexps)
    cals = fieldexps[fieldexps['flavor'] == flavor]
    notcal  = fieldexps[fieldexps['flavor'] != flavor]
    if len(cals) > 1:
        n_fixscreen_cal = len(set(cals[cals['FixedScreen']]['EXPOSURE']))
        if n_fixscreen_cal > 0:
            cals = cals[~cals['FixedScreen']]
            if len(cals) > 0:
                splog.info(f'Dropping {n_fixscreen_cal} {flavor} taken with fixed Dome flat screen')
            else:
                splog.info(f'Warning: No Windscreen {flavor}... Keeping {flavor} taken with fixed Dome flat screen')
    if (len(cals) > 0):
        fieldexps = vstack([cals,notcal])
    return fieldexps

def spplan_findrawdata(inputdir):
    """ Get the list of raw frames """
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

def get_FixedScreen_master_cal(flats, arcs):
    fflats = flats[flats['FixedScreen']].copy()
    farcs = arcs[arcs['FixedScreen']].copy()
    fflats.sort('TAI')
    farcs.sort('TAI')
    for arc in farcs:
        idx   = np.where(fflats['fieldid'].data == arc['fieldid'].data)[0]
        mflat = fflats[idx]
        if len(mflat) > 1:
            idx   = find_nearest(mflat['TAI'].data, arc['TAI'].data)
            mflat = mflat[idx]
            return mflat, arc
    return None, None

def get_FieldTaiMatch_master_cal(flats, arcs):
    marc  = arcs[0]
    idx   = np.where(flats['fieldid'].data == marc['fieldid'].data)[0]
    mflat = flats[idx]
    if len(mflat) > 1:
        idx   = find_nearest(mflat['TAI'].data, marc['TAI'].data)
        mflat = mflat[idx]
    elif len(mflat) == 0:
        ffields = np.unique(flats['fieldid'].data)
        afields = np.unique(arcs['fieldid'].data)
        mfields = np.intersect1d(ffields, afields)
        if len(mfields) > 0:
            mflat = flats[[fieldid in mfields for fieldid in flats['fieldid']]]
            marc  = arcs[[fieldid in mfields for fieldid in arcs['fieldid']]]
            mflat.sort('TAI')
            marc.sort('TAI')
            marc = marc[0]
            idx   = np.where(mflat['fieldid'].data == marc['fieldid'].data)[0]
            mflat = mflat[idx]
            if len(mflat) >1:
                idx   = find_nearest(mflat['TAI'].data, marc['TAI'].data)
                mflat = mflat[idx]
        else:
            idx   = find_nearest(flats['TAI'].data, marc['TAI'].data)
            mflat = flats[idx]
    return mflat, marc
    

def get_master_cal(allexps,dropMaster=True, obs='APO', mjd=''):
    allexps['flavor'] = allexps['flavor'].astype(object)
    if (obs.upper() == 'LCO'):
        allexps = mark_LCOfixedScreenCals(allexps)
    flats = allexps[np.where(allexps['flavor'].data == 'flat')[0]].copy()
    arcs  = allexps[np.where(allexps['flavor'].data == 'arc')[0]].copy()
    flats.sort('TAI')
    arcs.sort('TAI')
    if len(arcs) == 0 or len(flats) == 0:
        if len(arcs) == 0:
            splog.info(f'Warning: No Valid arcs (mjd:{mjd})')
        if len(flats) == 0:
            splog.info(f'Warning: No Valid Flats (mjd:{mjd})')
        return (None if not dropMaster else allexps)
        
    if (obs.upper() == 'LCO'):
        mflat, marc = get_FixedScreen_master_cal(flats, arcs)
        if mflat is None:
            mflat, marc = get_FieldTaiMatch_master_cal(flats, arcs)
    else:
        mflat, marc = get_FieldTaiMatch_master_cal(flats, arcs)

    if len(mflat) > 0:
        flats = allexps[np.where(allexps['EXPOSURE'].data == mflat['EXPOSURE'].data)[0]]
        arcs  = allexps[np.where(allexps['EXPOSURE'].data == marc['EXPOSURE'].data)[0]]
        flats['flavor']  = 'TRACEFLAT'
        arcs['flavor']   = 'TRACEARC'

    if dropMaster:
        drop = np.where(allexps['EXPOSURE'].data == flats['EXPOSURE'].data[0])[0]
        allexps.remove_rows(drop)
        drop = np.where(allexps['EXPOSURE'].data == arcs['EXPOSURE'].data[0])[0]
        allexps.remove_rows(drop)
    allexps = vstack([flats, arcs, allexps])

    allexps['flavor'] = allexps['flavor'].astype(str)
    return(allexps)

def build_exps(i, mj, mjdlist, OBS, rawdata_dir, ftype, spplan_Trace=False, no_remote=True,
                legacy=False, plates=False, fps=False, lco=False, release='sdsswork',
                verbose=True, no_dither=False):
    thismjd = int(mj)

    if OBS == 'APO':
        if thismjd  ==  59560:
            splog.info(f'Skipping {thismjd} for FPS Commissioning')
            return [],None ##FPS Commissioning
        if thismjd in [59760,59755,59746,59736,59727,59716,59713]:
            splog.info(f'Skipping {thismjd} for 6450Ang Feature')
            return [],None #6450 Feature:

    
    inputdir = ptt.join(rawdata_dir, mj)
    sphdrfix = Sphdrfix(mj, fps=ftype.fps, obs=OBS)
    splog.info('Data directory '+inputdir)

    if (len(mjdlist) == 1) and (not spplan_Trace):
        wait = check_transfer(OBS, mj)
        if int(mj) < 59148:
            wait = False
        i = 0
        while wait:
            if i <= 15:
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
            if spplan_Trace:
                if FLAVOR is None:
                    FLAVOR = ''
                if FLAVOR.lower() not in ['arc','flat','calibration']:
                    continue
            #-----------
            # Rename 'target' -> 'science', and 'calibration' -> 'arc'
            #try:
            if True:
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
#            except:
#                splog.info(f'Skipping '+ptt.basename(f)+' with failed fits header')
#                continue
                
            reject = Reject(f, hdr)
            if reject.check():
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
                              ' flavor '+FLAVOR+ ' inconsistent with map name for ' +
                              ptt.basename(f))
                    continue
                
                if len(MAPNAME) <=4:
                    """
                     MAPNAME should be of the form '000000-51683-01'.
                     If it only contains the FIELDID ;; (for MJD <= 51454),
                     then find the actual plug-map file.
                    """
                    plugfile = find_plPlugMapM('*', str(int(MAPNAME)).zfill(4), release=release, no_remote=no_remote)
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
                offset_ra = getcard(hdr,'OFFRA', default = 0, noNaN=True)
                offset_dec = getcard(hdr,'OFFDEC', default = 0, noNaN=True)
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
                    confile = get_confSummary(CONFNAME, obs=OBS, release=release,
                                              sort=False, filter=False, no_remote=no_remote)
                    last_conf = confile
                if len(confile) == 0:
                    splog.info(f'Warning: Invalid or Missing confSummary for {CONFNAME} for '+str(hdr['EXPOSURE']).strip())
                    if (fieldid != 0) and (fieldid != -999):
                        continue
                    
                if FLAVOR.lower() == 'science':
                    if 'is_dithered' in confile.meta.keys():
                        if (confile.meta['is_dithered'] == '1') or (confile.meta['parent_configuration'].strip() != '-999'):
                            dither = 'T'
                            if no_dither:
                                splog.info('Warning: Skipping Robot Dither Exposure '+str(hdr['EXPOSURE']).strip())
                                continue
                        
                    if (offset_ra !=0) or (offset_dec !=0):
                        expid = str(hdr['EXPOSURE']).strip()
                        if 59864 <= thismjd <= 60091 and OBS == 'APO' and offset_ra == 0 and offset_dec == -0.1:
                            pass
                            #splog.warning(f"Warning: Offset of dRA={round(offset_ra,2)} dDec={round(offset_dec,2)} on Science Frame {expid} due to GFA Offset correction")
                            # Offest implemented to correct for offset between GFAs and Wok
                        elif 60300 <= thismjd <= 60307 and OBS == 'APO' and 362988 <= int(expid) <= 363341:
                            # Telescope offset not reset after dithers...
                            splog.warning(f'Warning: Telescope offset of dRA={round(offset_ra,2)} dDec={round(offset_dec,2)} on Science Frame {expid}')
                        else:
                            dither = 'T'
                            if no_dither:
                                splog.info(f'Warning: Skipping Telescope Dither Exposure {expid} (dRA={round(offset_ra,2)}, dDec={round(offset_dec,2)})')
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
                        'alt':[float(getcard(hdr, 'ALT', default=90.0))],
                        'fieldid':[fid],
                        'confid':[str(CONFNAME)],
                        'CONFNAME':[str(CONFNAME)],
                        'shortname':[ptt.basename(f)]})

            allexps = vstack([allexps,exp])
    return allexps, ftype

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def pair_ccds(ftype, fieldexps, clobber=False, override_manual=False, OBS='LCO'):
    if ftype.legacy:
        shape = (4,)
        cams = ['b1', 'b2', 'r1', 'r2']
    elif OBS == 'LCO':
        shape = (2,)
        cams = ['b2','r2']
    elif OBS == 'APO':
        shape = (2,)
        cams = ['b1','r1']
    splog.info(f"Grouping Matching exposures for {','.join(cams)}")
    fieldexps.add_column(Column('', dtype=object, name='name', shape=shape))
    for i, row in enumerate(fieldexps):
        tnames = fieldexps[np.where(fieldexps['EXPOSURE'].data == row['EXPOSURE'])[0]]['shortname'].data
        tnames = np.char.replace(tnames, '.gz', '')
        tnames = list(set(tnames))
        names = []
        for cam in cams:
            match = False
            for tn in tnames:
                if cam in tn:
                    names.append(tn)
                    match = True
                    break
            if not match:
                splog.error(f"ERROR: Missing Frame: sdR-{cam}-{str(row['EXPOSURE']).zfill(8)}")
                names.append('')
        while len(names) < shape[0]:
            names.append('')
        fieldexps[i]['name'] = np.asarray(names)
    names = fieldexps['name'].data
    names = np.stack(names.tolist(), axis=0)
    names = names.astype(str)
    fieldexps.remove_column('name')
    fieldexps.add_column(Column(names, dtype=str, name='name', shape=shape))
    fieldexps = unique(fieldexps, keys=['EXPOSURE','flavor'])

    fieldexps = fieldexps['confid','fieldid','mjd','mapname','flavor','exptime','name']
    return fieldexps

def write_plan(planfile, fieldexps, meta={}, clobber=False, override_manual=False):
    fieldexps.meta=meta
    
    if ptt.exists(planfile):
        if clobber is False:
            splog.info('WARNING: Will not over-write plan file: ' + ptt.basename(planfile))
            return
        else:
            test = read_table_yanny(planfile, 'SPEXP')
            if not override_manual:
                if 'manual' in test.meta.keys():
                    if test.meta['manual'] == 'T':
                        splog.info('WARNING: Will not over-write manual plan file: ' + ptt.basename(planfile))
                        return
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
    return

def spplan2d():
    # topdir=None, run2d=None, mjd=None, mjdstart=None, mjdend=None,
    #          field= None, fieldstart = None, fieldend=None,
    #          matched_flats=False, nomatched_arcs=False, lco=False, minexp=1,
    #          clobber=False, release='sdsswork', logfile=None, no_remote=True,
    #          legacy=False, plates=False, fps=True, no_commissioning=False, no_dither=False,
    #          single_flat=False, single_arc=False, override_manual=False, manual_noarc=False,
    #          verbose = False, returnlist=False, **extra_kwds):
    
    logfile = config.pipe['plan.daily.dailyplan_logfile']
    if logfile is not None:
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    splog.info('spplan2d started at '+time.ctime())

    filt_field = config.pipe['fmjdselect.field']
    fieldstart = config.pipe['fmjdselect.fieldstart']
    fieldend = config.pipe['fmjdselect.fieldend']


    if config.pipe['fmjdselect.obs'].lower() == 'lco':
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
        OBS = 'LCO'
        if mjdstart is None:
            mjdstart = 60000
        elif int(mjdstart) <  60000:
            mjdstart = 60000
        mjdstart = int(mjdstart)
    else:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_N'
        OBS = 'APO'
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
                       start=config.pipe['fmjdselect.mjdstart'], 
                       end=config.pipe['fmjdselect.mjdend'])
    nmjd = len(mjdlist)
    splog.info(f'Number of MJDs = {nmjd}')
    if nmjd == 0:
        splog.info('No Valid MJDs')
        return None
    #---------------------------------------------------------------------------
    # Loop through each input MJD directory

    if returnlist:
        plans_list=[]
    dithered_pmjds = []

    legacy = config.pipe['SDSS_Generation.legacy']
    plates = config.pipe['SDSS_Generation.plates']
    fps = config.pipe['SDSS_Generation.fps']
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

            if traceflat:
                allexps = get_master_cal(allexps, dropMaster = False, obs= OBS, mjd=mj)
                if allexps is None:
                    continue
            for field in list(dict.fromkeys(allexps[fieldmap_col].data)):
                if ftype.legacy or ftype.plates:
                    ftype_exp = Fieldtype(mjd=mj)
                else:
                    ftype_exp = Fieldtype(fieldid=field_to_string(field), mjd=mj)
                if ftype.legacy is not ftype_exp.legacy:
                    splog.info(f'Warning: Skipping Legacy plate {field_to_string(field)} from non-Legacy MJD')
                    continue
                elif ftype.plates is not ftype_exp.plates:
                    splog.info(f'Warning: Skipping SDSS-V plate {field_to_string(field)} from non-SDSS-V Plate MJD')
                    continue
                elif ftype_exp.commissioning is True and no_commissioning is True:
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
                if not mjd_match(ftest, mjd=filt_field, mjdstart=fieldstart, mjdend=fieldend):
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
                elif nsci < minexp:
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

                if returnlist:
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
        # topdir=None, run2d=None, mjd=None, mjdstart=None, mjdend=None,
        #      field= None, fieldstart = None, fieldend=None, lco=False,
        #      clobber=False, logfile=None, override_manual=False,
        #      legacy=False, plates=False, plate_epoch = False,
        #      daily = False, plans=None, **extra_kwds):
    
    if not config.pipe['plan.daily.plate_epoch']:
        daily = True

    logfile = config.pipe['plan.daily.dailyplan_logfile']
    if config.pipe['plan.daily.skip2d']:
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

    if plans is not None:
        field = config.pipe['fmjdselect.field']
        if field is None:
            field = []
        field.extend([ptt.basename(x).split('-')[1] for x in np.atleast_1d(plans)])
    
    afc = Field(topdir, run2d, '*')
    fieldlist = get_dirs(ptt.dirname(afc.dir()), field = True,
                         match=field, start=fieldstart, end=fieldend)
    splog.info('Number of field directories = '+ str(len(fieldlist)))

    # Loop through each input configuration directory
    for fielddir in fieldlist:
        try: 
            fieldid = int(ptt.basename(fielddir))
        except:
            continue
        fc = Field(topdir, run2d, fielddir, mjd=mjd)
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
                test = mjd_match(row['mjd'], mjd=config.pipe['fmjdselect.mjd'], 
                                 mjdstart=config.pipe['fmjdselect.mjdstart'], 
                                 mjdend=config.pipe['fmjdselect.mjdend'])
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
                
                
