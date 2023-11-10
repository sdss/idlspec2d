#!/usr/bin/env python3

from splog import Splog
from os import getenv, makedirs, rename
import os.path as ptt
from glob import glob
import time
from field import field_to_string, Fieldtype
from astropy.table import Table, vstack, Column, unique
from astropy.io import fits
from collections import OrderedDict
from pydl.pydlutils.yanny import read_table_yanny, yanny, write_table_yanny
import subprocess
import numpy as np
import argparse
from load_module import load_module
from load_module import load_env
from spplan_epoch import get_dirs
from uubatchpbs import mjd_match
from GetconfSummary import find_confSummary, find_plPlugMapM, get_confSummary

from sdss_access.path import Path
from sdss_access import Access

from pydl import __version__ as pydlVersion
from sdss_access import __version__ as saver
from tree import __version__ as treever


splog = Splog()

SDSSCOREVersion = getenv('SDSSCORE_VER', default= '')
speclogVersion = subprocess.getoutput("speclog_version")
idlspec2dVersion = subprocess.getoutput("idlspec2d_version")
idlutilsVersion = subprocess.getoutput("idlutils_version")


def getcard(hdr, card, default=None, noNaN=False):
    try:
        if hdr.count(card) > 0:
            if type(hdr[card]) is str:
                hdr[card] = hdr[card].strip().replace("'","")
                if noNaN is True:
                    if hdr[card].strip().upper() == 'NAN':
                        hdr[card] = default
                return(hdr[card])
            else:
                return(hdr[card])
        else:
            return(default)
    except:
        return(default)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def get_alt_cal(fieldexps, allexps, flav='arc', single_cal=False):
    cals  = allexps[np.where(allexps['flavor'].data == flav)[0]].copy()
    idx_n0 = np.where(cals['fieldid'].data != field_to_string(0))[0]

    if len(idx_n0) != 0:
        cals = cals[idx_n0]

    if len(cals) == 0:
        return(fieldexps)



    if single_cal is True:
        texp = np.nanmean(fieldexps['TAI'].data)
        idx = find_nearest(cals['TAI'].data, texp)
        cal = cals[idx]
        idx = np.where(cals['EXPOSURE'].data == cal['EXPOSURE'].data)[0]
        cals = cals[idx]
    cals['fieldid'] = fieldexps[0]['fieldid'].data
    fieldexps = vstack([cals,fieldexps])

    return(fieldexps)

def get_master_cal(allexps):
    allexps['flavor'] = allexps['flavor'].astype(object)
    flats = allexps[np.where(allexps['flavor'].data == 'flat')[0]].copy()
    arcs  = allexps[np.where(allexps['flavor'].data == 'arc')[0]].copy()
    marc  = arcs[0]
    idx   = np.where(flats['fieldid'].data == marc['fieldid'].data)[0]
    mflat = flats[idx]
    if len(flats) > 1:
        idx   = find_nearest(mflat['TAI'].data, marc['TAI'].data)
        mflat = mflat[idx]
    if len(flats) > 0:
        flats = allexps[np.where(allexps['EXPOSURE'].data == mflat['EXPOSURE'].data)[0]]
        arcs  = allexps[np.where(allexps['EXPOSURE'].data == marc['EXPOSURE'].data)[0]]
        flats['flavor']  = 'TRACEFLAT'
        arcs['flavor']   = 'TRACEARC'

    drop = np.where(allexps['EXPOSURE'].data == flats['EXPOSURE'].data[0])[0]
    allexps.remove_rows(drop)
    drop = np.where(allexps['EXPOSURE'].data == arcs['EXPOSURE'].data[0])[0]
    allexps.remove_rows(drop)
    allexps = vstack([flats, arcs, allexps])

    allexps['flavor'] = allexps['flavor'].astype(str)
    return(allexps)


class Sphdrfix:
    def __init__(self, mjd, fps=False, obs='APO', release=None, no_remote=True):
        self.mjd = str(mjd)
        self.sphdrfix_table = None
        self.fps = fps
        self.obs = obs.lower()
        self.release = release
        self.no_remote = no_remote
        self.read_sphdrfix()
        
    def read_sphdrfix(self):
        if self.release is not None:
            path = Path(release=self.release, preserve_envvars=True)
            path_options = {'mjd':self.mjd}
            if path.exists('sdHdrFix', **path_options):
                reportfile = path.full('sdHdrFix', **path_options)
            elif path.exists('sdHdrFix', **path_options, remote=(not self.no_remote)):
                access = Access(release=self.release)
                reportfile = path.full('sdHdrFix', **path_options)
                access.remote()
                access.add('sdHdrFix', **path_ops)
                access.set_stream()
                valid = access.commit()
                if valid is False:
                    return
            else:
                return
        elif not self.fps:
            speclog_dir = getenv('SPECLOG_DIR')
            if speclog_dir is None:
                splog.info('ERROR: Must set environment variabel SPECLOG_DIR')
                exit(1)
            reportfile = ptt.join(speclog_dir, self.mjd, 'sdHdrFix-'+self.mjd+'.par')
        else:
            speclog_dir = getenv('SDHDRFIX_DIR')
            if speclog_dir is None:
                splog.info('ERROR: Must set environment variabel SDHDRFIX_DIR')
                exit(1)
            reportfile = ptt.join(speclog_dir, self.obs, 'sdHdrfix','sdHdrFix-'+self.mjd+'.par')
        if ptt.exists(reportfile):
            self.sphdrfix_table = read_table_yanny(reportfile, 'OPHDRFIX')
            self.sphdrfix_table.convert_bytestring_to_unicode()

    def fix(self, infile, hdr):
        fileroot = ptt.basename(infile).split('.')[0]
        wfileroot = fileroot.split('-')
        wfileroot = '-'.join([wfileroot[0], '??', wfileroot[-1]])
        if self.sphdrfix_table is not None:
            for row in self.sphdrfix_table:
                if (row['fileroot'] == fileroot) or (row['fileroot'] == wfileroot):
                    hdr[row['keyword']] = row['value']
        if getcard(hdr,'QUALITY') is None:
            hdr['QUALITY'] = 'excellent'
        return hdr
    
def get_key(fp):
    filename = ptt.splitext(ptt.splitext(ptt.basename(fp))[0])[0]
    int_part = filename.split('-')[2]
    try:
        return int(int_part)
    except:
        return int(ptt.splitext(int_part)[0])
    

#from tqdm import tqdm
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
        
        
def spplanTrace(topdir=None, run2d=None, mjd=None, mjdstart=None, mjdend=None,
             lco=False, clobber=False, release='sdsswork', logfile=None, no_remote=True,
             legacy=False, plates=False, override_manual=False, 
             verbose = False, **extra_kwds):
    
    if logfile is not None:
        splog.open(logfile=logfile, logprint=False)
        splog.info('Log file '+logfile+' opened '+ time.ctime())
    else:
        if 'splog' in extra_kwds.keys():
            splog = extra_kwds['splog']
        else:
            splog = globals()['splog']
        splog.info('spPlanTarget started at '+time.ctime())

    if lco:
        BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
        OBS = 'LCO'
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
    
   #----------
   # Create a list of the MJD directories (as strings)
    mjdlist = get_dirs(rawdata_dir, subdir='', pattern='*', match=mjd, start=mjdstart, end=mjdend)
    nmjd = len(mjdlist)
    splog.info(f'Number of MJDs = {nmjd}')
    
    plateflavors = ['BHM', 'BHM&MWM', 'EBOSS', 'BOSS']
    #---------------------------------------------------------------------------
    # Loop through each input MJD directory

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
        sphdrfix = Sphdrfix(mj, fps=ftype.fps, obs=OBS)

#        if ftype.legacy or ftype.plates:
#            plugdir = ptt.join(speclog_dir, mj)
        splog.info('----------------------------')
        splog.info(f'MJD: {mj} {ftype} ({i+1} of {len(mjdlist)})')
        splog.info('Data directory '+inputdir)


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
                if FLAVOR.lower() not in ['arc','flat','calibration']:
                    continue
                #-----------
                # Rename 'target' -> 'science', and 'calibration' -> 'arc'
    
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
                        splog.infof(f'Skipping {FLAVOR} file '+ptt.basename(f)+' with closed Hartmann Doors')
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
            splog.info('No Calibration Frames for '+mj)
        else:
            if ftype.legacy or ftype.plates:
                fieldmap_col = 'mapname'
            else:
                fieldmap_col = 'fieldid'

            
            manual = 'F'
            allexps = get_master_cal(allexps)
            planfile = 'spPlanTrace-' + mj + '_'+OBS+'.par'
            planfile = ptt.join(topdir, run2d, 'trace', str(thismjd), planfile)
                
                
            if ftype.legacy:
                shape = (4,)
            else:
                shape = (2,)
                
            allexps.add_column(Column('', dtype=object, name='name', shape=shape))

            for i, row in enumerate(allexps):
                tnames = allexps[np.where(allexps['EXPOSURE'].data == row['EXPOSURE'])[0]]['shortname'].data
                tnames = np.char.replace(tnames, '.gz', '')
                names = []
                for cam in ['b1', 'b2', 'r1', 'r2']:
                    for tn in tnames:
                        if cam in tn:
                            names.append(tn)
                while len(names) < shape[0]:
                    names.append('')
                allexps[i]['name'] = np.asarray(names)
            names = allexps['name'].data
            names = np.stack(names.tolist(), axis=0)
            names = names.astype(str)
            allexps.remove_column('name')
            allexps.add_column(Column(names, dtype=str, name='name', shape=shape))
            allexps = unique(allexps, keys='EXPOSURE')

            allexps = allexps['confid','fieldid','mjd','mapname','flavor','exptime','name']
        
                
            allexps.meta=OrderedDict({
                            'MJD':              mj                       +"   # Modified Julian Date",
                            'OBS':              OBS                      +"   # Observatory",
                            'RUN2D':            run2d                    +"   # 2D reduction name",
                            'idlspec2dVersion': "'"+idlspec2dVersion+"'" +"   # idlspec2d Version when building plan",
                            'idlutilsVersion':  "'"+idlutilsVersion+"'"  +"   # idlutils Version when building plan",
                            'pydlVersion':      "'"+pydlVersion+"'"      +"   # Version of pydl when building plan",
                            'speclogVersion':   "'"+speclogVersion+"'"   +"   # speclog Version when building plan",
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
            allexps.convert_unicode_to_bytestring()
            write_table_yanny(allexps, planfile,tablename='SPEXP', overwrite=clobber)
        del allexps
    splog.info('----------------------------')
    splog.info('Successful completion of spplanTrace at '+ time.ctime())

    if logfile is not None:
        splog.close()
    return
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Produces spPlanTrace')
    

    General = parser.add_argument_group(title='General', description='General Setup Options')
    General.add_argument('--module',         help='Module file to load for run')
    General.add_argument('--topdir',         help='')
    General.add_argument('--run2d',          help='Run2d to override module or environmental variable')
    General.add_argument('--lco',            help='Build Run files for LCO', action='store_true')
    General.add_argument('--logfile',        help='Optional logfile (Including path)')
    General.add_argument('--verbose',        help='Provide information about nonutlized frames')
    General.add_argument('-c', '--clobber',  help='overwrites previous plan file', action='store_true')
    General.add_argument('--release',        help='sdss_access data release (defaults to sdsswork), required if you do not have proprietary access, otherwise see https://sdss-access.readthedocs.io/en/latest/auth.html#auth', default='sdsswork')
    General.add_argument('--remote',         help='allow for remote access to data using sdss-access', action='store_true')
    General.add_argument('--override_manual',help='Override/clobber manually edited plan', action='store_true')

    
    
    Filter_grp = parser.add_argument_group(title='MJD/Field Filtering', description='MJD/Field Filtering Options')
    Filter_grp.add_argument('--mjd',  nargs='*',  help = 'Use data from these MJDs.')
    Filter_grp.add_argument('--mjdstart',         help = 'Starting MJD')
    Filter_grp.add_argument('--mjdend',           help = 'Ending MJD')

    Filter_grp.add_argument('--legacy',           help = 'Include legacy (BOSS/eBOSS) plates',   action='store_true')
    Filter_grp.add_argument('--plates',           help = 'Include SDSS-V plates',                action='store_true')
    Filter_grp.add_argument('--fps',              help = 'Include FPS Fields',                   action='store_true')
    Filter_grp.add_argument('--sdssv',            help = 'Include both SDSS-V Fields & Plates',  action='store_true')
    
    
    args = parser.parse_args()

    if args.sdssv:
        args.fps = True
        args.plates = True

    if args.release != 'sdsswork':
        if args.release not in Access().get_available_releases():
            parser.exit(status=0, message='ERORR: '+args.release+' is not a valid release')
    else:
        if args.remote is True:
            try:
                Access().remote()
            except:
                parser.exit(status=0, message='ERROR: No netrc file found. see https://sdss-access.readthedocs.io/en/latest/auth.html#auth')

    if args.module is not None:
        module = load_module()
        module('purge')
        module('load', args.module)
        if args.run2d is None:
            args.run2d = load_env('RUN2D')
        if args.topdir is None:
            args.topdir = load_env('BOSS_SPECTRO_REDUX')

    args.no_remote = not args.remote

    idlspec2dVersion = subprocess.run(['idlspec2d_version'], stdout=subprocess.PIPE).stdout.decode('utf-8')
    idlutilsVersion  = subprocess.run(['idlutils_version'],  stdout=subprocess.PIPE).stdout.decode('utf-8')
    speclogVersion   = subprocess.run(['speclog_version'],   stdout=subprocess.PIPE).stdout.decode('utf-8')
    SDSSCOREVersion  = load_env('SDSSCORE_VER')


    idlspec2dVersion = idlspec2dVersion.replace('\n', '')
    idlutilsVersion = idlutilsVersion.replace('\n', '')
    speclogVersion = speclogVersion.replace('\n', '')
    spplanTrace(**vars(args))
