#!/usr/bin/env python3

import argparse
import sys
import os.path as ptt
from os import getenv, makedirs, remove
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, vstack, join, unique, setdiff
from fieldlist import wwhere, plot_sky, fieldlist, get_lastline
from healpy import ang2pix
from pydl.pydlutils.spheregroup import spheregroup
import time
from datetime import timedelta
import warnings
from field import field_to_string
from merge_dm import merge_dm
from glob import glob
from splog import Splog
import gc
splog = Splog()
run2d_warn = True

from sdss_semaphore.targeting import TargetingFlags


def read_zans(sp1d_dir, field, mjd):
    zansfile = ptt.join(sp1d_dir,'spZbest-' + field + '-' + mjd + '.fits')
    if ptt.exists(zansfile):
        splog.log('Reading Zbest file: '+ptt.basename(zansfile))
        try:
            zans = Table(fits.getdata(zansfile))
            for key in zans.colnames:
                zans.rename_column(key,key.upper())
        except:
            splog.log('Error reading '+zansfile)
            zans = None
    else:
        splog.log(zansfile +' missing')
        zans = None
    return(zans)

def read_zline(sp1d_dir, field, mjd):
    zlinefile = ptt.join(sp1d_dir,'spZline-' + field + '-' + mjd + '.fits')
    if ptt.exists(zlinefile):
        splog.log('Reading Zline file: '+ptt.basename(zlinefile))
        try:
            zline = Table(fits.getdata(zlinefile))
            for key in zline.colnames:
                zline.rename_column(key,key.upper())
        except:
            splog.log('Error reading '+zlinefile)
            zline = None
    else:
        splog.log(zlinefile +' missing')
        zline = None
    return(zline)

def read_xcsao(sp1d_dir, field, mjd, spAll):
    XCSAOfile = ptt.join(sp1d_dir,'spXCSAO-' + field + '-' + str(mjd) + '.fits')
    if ptt.exists(XCSAOfile):
        splog.log('Reading XCSAO file: '+ptt.basename(XCSAOfile))
        try:
            XCSAO_tab = Table(fits.getdata(XCSAOfile,1))
            for key in XCSAO_tab.colnames:
                XCSAO_tab.rename_column(key,key.upper())
            XCSAO_tab.remove_columns(['EBV','FIBERID_LIST','OBJID','FIELDID'])
            
            renames = {'RV':'XCSAO_RV','ERV':'XCSAO_ERV','R':'XCSAO_RXC','TEFF':'XCSAO_TEFF',
                       'ETEFF':'XCSAO_ETEFF', 'LOGG':'XCSAO_LOGG', 'ELOGG':'XCSAO_ELOGG',
                       'FEH':'XCSAO_FEH', 'EFEH':'XCSAO_EFEH'}
            
            for col in renames.keys():
                XCSAO_tab.rename_column(col,renames[col])
            for col in XCSAO_tab.colnames:
                if col == 'TARGET_INDEX':
                    continue
                if col in spAll.colnames:
                    XCSAO_tab.remove_column(col)
            if 'TARGET_INDEX' not in XCSAO_tab.colnames:
                splog.log('Error reading '+XCSAOfile)
                XCSAO_tab = None
        except:
            splog.log('Error reading '+XCSAOfile)
            XCSAO_tab = None
    else:
        splog.log(XCSAOfile +' missing')
        XCSAO_tab = None
    return(XCSAO_tab)

def read_fibermap(field_dir, field, mjd, epoch=False, allsky=False):
    if epoch:
        spfieldfile = ptt.join(field_dir, 'epoch', 'spField-'+field+'-'+mjd+'.fits')
    else:
        spfieldfile = ptt.join(field_dir, 'spField-'+field+'-'+mjd+'.fits')
    if allsky is True:
        spfieldfile = spfieldfile.replace('spField-','spFullsky-')
    if ptt.exists(spfieldfile):
        splog.log('Reading Fibermap: '+ptt.basename(spfieldfile))
        try:
            fibermap = Table(fits.getdata(spfieldfile, 5))
            hdr = fits.getheader(spfieldfile,0)
            fibermap['SPEC_FILE'] = np.char.add(np.char.add('spec-'+field+'-'+mjd+'-', np.char.strip(np.asarray(fibermap['CATALOGID']).astype(str))),'.fits')
            renames = {'XFOCAL_LIST':'XFOCAL','YFOCAL_LIST':'YFOCAL','CARTON_TO_TARGET_PK_LIST':'CARTON_TO_TARGET_PK','ASSIGNED_LIST':'ASSIGNED',
                       'ON_TARGET_LIST':'ON_TARGET','VALID_LIST':'VALID','DECOLLIDED_LIST':'DECOLLIDED','ICATALOGID':'CATALOGID',
                       'GAIA_ID_DR2': 'GAIA_ID', 'MJDLIST': 'MJD_LIST','PROGRAM':'PROGRAMNAME'}
        
            for key in fibermap.colnames:
                fibermap.rename_column(key,key.upper())
            for key in renames.keys():
                if key in fibermap.columns:
                    if renames[key] in fibermap.columns:
                        fibermap.remove_column(renames[key])
                    fibermap.rename_column(key,renames[key])

            ra  = np.asarray(fibermap['RACAT']).astype(float)
            dec = np.asarray(fibermap['DECCAT']).astype(float)
           
            idx = np.where((np.isfinite(ra)) & (np.isfinite(dec)))[0]
            hp = np.zeros(len(ra), dtype = int)
            hp[idx] = ang2pix(128, ra[idx], dec[idx], lonlat=True)
            fibermap.add_column(hp, name = 'HEALPIX')
            fibermap.add_column(np.floor(fibermap['HEALPIX']/1000).astype(int), name = 'HEALPIXGRP')


            fibermap.add_column(-999, name = 'CATALOGID_V0')
            fibermap.add_column(-999, name = 'CATALOGID_V0P5')
            catid = fibermap['CATALOGID']
            cv0   = fibermap['CATALOGID_V0']
            cv0p5 = fibermap['CATALOGID_V0P5']
            iv0   = np.where(wwhere(fibermap['CATVERSION'].data,'0.0*'))[0]
            iv0p1   = np.where(wwhere(fibermap['CATVERSION'].data,'0.1*'))[0]
            iv0p5 = np.where(wwhere(fibermap['CATVERSION'].data,'0.5*'))[0]
            cv0[iv0]     = catid[iv0]
            cv0[iv0p1]   = catid[iv0p1]
            cv0p5[iv0p5] = catid[iv0p5]

            HEALPIX_PATH = np.char.add(np.char.add('$MWM_HEALPIX/', np.asarray(fibermap['HEALPIXGRP']).astype(str)),'/')
            HEALPIX_PATH = np.char.add(np.char.add(HEALPIX_PATH, np.asarray(fibermap['HEALPIX']).astype(str)),'/boss/')
            HEALPIX_PATH = np.char.add(HEALPIX_PATH, fibermap['SPEC_FILE'])
            idx = np.where(catid == 0)[0]
            HEALPIX_PATH[idx] = ''
            fibermap.add_column(HEALPIX_PATH, name = 'HEALPIX_PATH')
        except:
            splog.log('Error reading '+spfieldfile)
            fibermap = Table()
    else:
        splog.log(spfieldfile +' missing')
        fibermap = Table()
    return(fibermap)

def read_spline(splinefile, skip_line=False, clobber=False):
    if clobber:
        exists = False
    else:
        exists = ptt.exists(splinefile)
        if not exists:
            splinefile = splinefile.replace('.gz','')
            exists = ptt.exists(splinefile)
        if not exists:
            splinefile = splinefile+'.gz'
            exists = ptt.exists(splinefile)
    if exists and not skip_line:
        try:
            splog.log('Reading spline file: '+ptt.basename(splinefile))
            spline = Table(fits.getdata(splinefile))
            if len(spline) == 0:
                remove(splinefile)
                spline = None
        except:
            splog.log('Failure opening '+ splinefile)
            spline = None
    else:
        spline = None
    return(spline)

def read_spall(spAllfile, clobber=False):
    if clobber:
        exists = False
    else:
        exists = ptt.exists(spAllfile)
        if not exists:
            spAllfile = spAllfile.replace('.gz','')
            exists = ptt.exists(spAllfile)
        if not exists:
            spAllfile = spAllfile+'.gz'
            exists = ptt.exists(spAllfile)
    if exists:
        try:
            splog.log('Reading spAll file: '+ptt.basename(spAllfile))
            spAll = Table(fits.getdata(spAllfile))
            if len(spAll) == 0:
                remove(spAllfile)
                spAll = None
        except:
            splog.log('Failure opening '+ spAllfile)
            spAll = None
    else:
        spAll = None
    return(spAll)




def oneField(row, field, mjd, skip_line=False, include_bad=False, legacy=False, dev=False,
            skip_specprimary=False, XCSAO=False, indir='.', clobber=False, epoch=False,
            merge_only=False, custom = None, allsky=None):
    field=field_to_string(field)
    mjd = str(mjd)
    row['RUN1D'] = row['RUN1D'].strip()
    if len(row['RUN1D']) == 0:
        row['RUN1D'] = row['RUN2D']
    if custom is None:
        if epoch is True:
            sp1d_dir =  ptt.join(indir, row['RUN2D'], field, 'epoch', row['RUN1D'])
            field_dir =  ptt.join(indir,row['RUN2D'], field)
        else:
            sp1d_dir =  ptt.join(indir, row['RUN2D'], field, row['RUN1D'])
            field_dir =  ptt.join(indir, row['RUN2D'],field)
    else:
        if epoch is True:
            sp1d_dir =  ptt.join(indir, row['RUN2D'], custom, 'epoch', row['RUN1D'])
            field_dir =  ptt.join(indir,row['RUN2D'], custom)
        else:
            sp1d_dir =  ptt.join(indir, row['RUN2D'], custom, row['RUN1D'])
            field_dir =  ptt.join(indir, row['RUN2D'],custom)
        field = custom
    dev = False
    spAllfile, spAlllitefile, splinefile, spAlldatfile = build_fname(indir, row['RUN2D'],
                                                                     field=field, mjd=mjd,
                                                                     dev=dev, epoch=epoch,
                                                                     custom=custom, allsky=allsky)
    
    spAll  = read_spall(spAllfile, clobber=clobber)
    spline = read_spline(splinefile, skip_line = skip_line, clobber=clobber)

    if merge_only:
        return({'spall':spAll, 'spline': spline})


    ############################
    #  Build SpAll and spAllline
    ############################
    zans = None
    if spAll is None:
        zans = read_zans(sp1d_dir, field, mjd)
        spAll = read_fibermap(field_dir, field, mjd, epoch=epoch, allsky=allsky)
        if len(spAll) == 0:
            spAll = None
            
        else:
            fieldlist_keys = {'CHUNK':'CHUNK', 'FIELDQUALITY':'FIELDQUALITY',
                              'FIELDSN2':'FIELDSN2', 'OBSERVATORY':'OBS',
                              'DEREDSN2':'DEREDSN2', 'DESIGNID':'DESIGNID',
                              'FIELD':'FIELD'}
            if 'FIELD' in spAll.colnames:
                spAll.remove_column('FIELD')
            
            if len(spAll) > 0:
                for key in fieldlist_keys.keys():
                    if allsky:
                        if key == 'OBSERVATORY': continue
                    if key in row.colnames:
                        spAll.add_column(row[key],name=fieldlist_keys[key])

            if XCSAO:
                XCSAO_tab = read_xcsao(sp1d_dir, field, mjd, spAll)
                if XCSAO_tab is not None:
                    spAll = join(spAll,XCSAO_tab,keys='TARGET_INDEX', join_type='left')

            if zans is None:
                return ({'spall':None, 'spline': None})
            for col in zans.colnames:
                if col == 'TARGET_INDEX':
                    continue
                if col in spAll.colnames:
                    zans.remove_column(col)
            spAll = join(spAll,zans,keys='TARGET_INDEX', join_type='left')

        if spAll is not None:
            if len(spAll) > 0:
                for key in spAll.colnames:
                    spAll.rename_column(key,key.upper())
                spAll['FIBER_RA']  = spAll['FIBER_RA'].data.astype(float)
                spAll['FIBER_DEC'] = spAll['FIBER_DEC'].data.astype(float)
                spAll = build_specobjid(spAll)
                
    if spline is None and skip_line is False:
        spline = read_zline(sp1d_dir, field, mjd)

    if spline is not None:
        if len(spline) > 0:
            for key in spline.colnames:
                spline.rename_column(key,key.upper())
    if spAll is not None:
        if len(spAll) == 0:
            spAll = None
    if spline is not None:
        if len(spline) == 0:
            spline = None
    return({'spall':spAll, 'spline': spline})



def build_specobjid(spAll):
    global run2d_warn

    lsh = lambda x,s: np.uint64(x)*np.uint64(2**s)


    run2ds = spAll['RUN2D'].data
    
    run2ds = np.char.strip(run2ds)
    rerun = np.zeros(len(run2ds))
    
    if 'SPECOBJID' not in spAll.colnames:
        spAll.add_column(np.uint64(0),name='SPECOBJID')
    else:
        spAll.remove_column('SPECOBJID')
        spAll.add_column(np.uint64(0),name='SPECOBJID')
    
    for row in spAll:
        try:
            if isinstance(row['RUN2D'].strip(),str):
                n,m,p = row['RUN2D'].strip().split('_')
                n = int(n[1:])
                m = int(m)
                p = int(p)
                run2d = (n-5)*10000 + m*100 + p
            elif row['RUN2D'].strip().isnumeric():
                run2d = int(run2d)
            else:
                if run2d_warn:
                    splog.log(f"WARNING: Unable to parse RERUN from {np.unique(run2ds)} for CAS-style SPECOBJID; Using 0 instead")
                    run2d_warn = False
                run2d = 0
        except:
            if run2d_warn:
                splog.log(f"WARNING: Unable to parse RERUN from {np.unique(run2ds)} for CAS-style SPECOBJID; Using 0 instead")
                run2d_warn = False
            run2d = 0
        specobjid = np.uint64(0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            specobjid |= lsh(row['FIELD'],50) | lsh(row['TARGET_INDEX'],38) | lsh(row['MJD']-50000,24) | lsh(run2d,10)
        row['SPECOBJID'] = specobjid

    return(spAll)


def specPrimary(spAll):
    t2 = time.time()
    # Determine the score for each object
    # 1) Prefer observations with positive SN_MEDIAN in r-band
    # 2) Prefer fieldQUALITY='good' over any other field quality
    # 3) Prefer observations with ZWARNING=0
    # 4) Prefer objects with larger SN_MEDIAN in r-band
    # ASBjuly2011: test against ZWARNING_NOQSO for GALAXY targets:
    splog.log(f'Assigning primaries (start: {time.ctime()})')
    zw_primtest = spAll['ZWARNING']

    wh_galtarget = np.where(wwhere(spAll['OBJTYPE'], 'GALAXY*'))[0]
    if len(wh_galtarget) > 0:
        zw_primtest[wh_galtarget] = spAll['ZWARNING_NOQSO'][wh_galtarget]
    
    if len(spAll['SN_MEDIAN'][0]) == 1:
        jfilt = 0 
    else:
        jfilt = 2

    score = (4 * (spAll['SN_MEDIAN'][:,jfilt] > 0) + 2*(wwhere(spAll['FIELDQUALITY'],'good*')) 
            + 1 * (zw_primtest == 0) + (spAll['SN_MEDIAN'][:,jfilt] > 0)) / max(spAll['SN_MEDIAN'][:,jfilt]+1.)

    ra  = np.asarray(spAll['RACAT']).astype(float)
    dec = np.asarray(spAll['DECCAT']).astype(float)

    idx = np.where((np.isfinite(ra)) & (np.isfinite(dec)))[0]
    ingroup, multgroup, firstgroup, nextgroup = spheregroup(ra[idx], dec[idx], 2.0 / 3600., chunksize=None)
    # Set the unique object IDs
    spAll[idx]['BOSS_SPECOBJ_ID'] = ingroup + 1

    for col in ['SPECPRIMARY','SPECBOSS', 'NSPECOBS']:
        if col not in spAll.colnames:
            spAll.add_column(0,name=col)
        else:
            spAll[col] = 0

    sprim = spAll['SPECPRIMARY']
    nspec = spAll['NSPECOBS']
    sboss = spAll['SPECBOSS']
    spAll['NSPECOBS'] = 1
    for j, fg in enumerate(firstgroup):
        if fg !=1:
            if multgroup[j] == 0:
                continue
            elif multgroup[j] == 1:
                sprim[idx[fg]] = 1
                nspec[idx[fg]] = 1
                sboss[idx[fg]] = 1
            else:
                indx = np.zeros(multgroup[j], dtype=int)
                indx[0] = fg
                for k in range(multgroup[j] -1):
                    indx[k+1] = nextgroup[indx[k]]
                ibest = np.argmax(score[idx][indx])
                sprim[idx[indx[ibest]]] = 1
                nspec[idx[indx[ibest]]] = multgroup[j]
                sboss[idx[indx[ibest]]] = 1
    splog.log('Time to assign primaries = '+str(timedelta(seconds = time.time() - t2)))
    return(spAll)


def specPrimary_sdssid(spAll):
    t2 = time.time()
    # Determine the score for each object
    # 1) Prefer observations with positive SN_MEDIAN in r-band
    # 2) Prefer fieldQUALITY='good' over any other field quality
    # 3) Prefer observations with ZWARNING=0
    # 4) Prefer objects with larger SN_MEDIAN in r-band
    # ASBjuly2011: test against ZWARNING_NOQSO for GALAXY targets:
    splog.log(f'Assigning primaries (start: {time.ctime()})')
    zw_primtest = spAll['ZWARNING']

    wh_galtarget = np.where(wwhere(spAll['OBJTYPE'], 'GALAXY*'))[0]
    if len(wh_galtarget) > 0:
        zw_primtest[wh_galtarget] = spAll['ZWARNING_NOQSO'][wh_galtarget]
    
    if len(spAll['SN_MEDIAN'][0]) == 1:
        jfilt = 0
    else:
        jfilt = 2

    score = (4 * (spAll['SN_MEDIAN'][:,jfilt] > 0) + 2*(wwhere(spAll['FIELDQUALITY'],'good*'))
            + 1 * (zw_primtest == 0) + (spAll['SN_MEDIAN'][:,jfilt] > 0)) / max(spAll['SN_MEDIAN'][:,jfilt]+1.)
    spAll.add_column(score,'SCORE')
    sdssids = np.unique(spAll['SDSS_ID'].values)
    
    spAll.add_column(0,name='SPECPRIMARY')
    spAll.add_column(0,name='SPECBOSS')
    spAll.add_column(0,name='NSPECOBS')
    for id in sdssids:
        idx = np.where(SpAll['SDSS_ID'].values == id)[0]
        SpAll[idx]['NSPECOBS'] = len(idx)
        primary = np.argmax(SpAll[idx]['SCORE'].values)
        spAll[idx[primary]]['SPECPRIMARY'] = 1
        spAll[idx[primary]]['SPECBOSS'] = 1
    splog.log('Time to assign primaries = '+str(timedelta(seconds = time.time() - t2)))
    return(spAll)


def build_custom_fieldlist(indir, custom, run2d, run1d):
    flist = Table()
    for spfield in glob(ptt.join(indir, run2d,custom,f'spFullsky-{custom}-*.fits')):
        hdr = fits.getheader(spfield,0)
        
        spdiagcomblog = ptt.join(indir, run2d,custom,f"spDiagcomb-{custom}-{str(hdr['RUNMJD'])}.log")
        if ptt.exists(spdiagcomblog):
            lastline = get_lastline(spdiagcomblog)
            if 'Successful completion' in lastline:
                #Case where this 1D log file completed, which is not a case that should ever occur
                STATUSCOMBINE = 'Done'
            else:
                #Case where this 1D log file isn't completed
                STATUSCOMBINE = 'RUNNING'
        else:
            STATUSCOMBINE = 'Pending'
        spDiag1dlog = ptt.join(indir, run2d, custom, run1d, f"spDiag1d-{custom}-{str(hdr['MJD'])}.log")
        if ptt.exists(spDiag1dlog):
            lastline = get_lastline(spDiag1dlog)
            if 'Successful completion' in lastline:
                #Case where this 1D log file completed, which is not a case that should ever occur
                STATUS1D = 'Done'
            else:
                #Case where this 1D log file isn't completed
                STATUS1D = 'RUNNING'
        else:
            STATUS1D = 'Pending'

        try:
            sn2_g1 = hdr['SPEC1_G']
        except:
            sn2_g1 = np.NaN
        try:
            sn2_i1 = hdr['SPEC1_I']
        except:
            sn2_i1 = np.NaN

        try:
            sn2_g2 = hdr['SPEC2_G']
        except:
            sn2_g2 = np.NaN
        try:
            sn2_i2 = hdr['SPEC2_I']
        except:
            sn2_i2 = np.NaN
        flist_row = Table({'RUN2D':[run2d], 'RUN1D':[run1d],
                           'PROGRAM':[custom], 'FIELD':[0], 'MJD': [hdr['MJD']],
                           'STATUS2D':['Done'], 'STATUSCOMBINE':[STATUSCOMBINE],
                           'STATUS1D':[STATUS1D],'FIELDQUALITY':['good'],
                           'FIELDSN2':[np.nanmin( np.array([sn2_g1,sn2_i1,sn2_g2,sn2_i2]))],
                           'OBS':['']})
    
        flist =  vstack([flist, flist_row])
    print(flist)
    return(flist)


def fieldmerge(run2d=getenv('RUN2D'), indir= getenv('BOSS_SPECTRO_REDUX'),
               skip_line=False, include_bad=False, legacy=False, skip_specprimary=False,
               lite=False, XCSAO=False, field=None, mjd=None, programs=None, clobber=False, dev=False,
               datamodel=None, line_datamodel=None, verbose=False, epoch =False, outroot=None,
               logfile=None, remerge_fmjd=None, merge_only=False, limit=None,
               custom=None, allsky=False, run1d=None):
               
    try:
        SDSSC2BV = TargetingFlags.meta['SDSSC2BV']
    except:
        SDSSC2BV = '1'
    if field is not None:
        field = field_to_string(field)
    if mjd is not None:
        mjd = str(mjd)
        
    if datamodel is None:
        datamodel = ptt.join(getenv('IDLSPEC2D_DIR'), 'datamodel', 'spall_dm.par')
    if line_datamodel is None:
        line_datamodel = ptt.join(getenv('IDLSPEC2D_DIR'), 'datamodel', 'spzline_dm.par')
 
    if epoch is True:
        fieldlist_file = ptt.join(indir, run2d, 'epoch', 'fieldlist-'+run2d+'.fits')
    else:
        fieldlist_file = ptt.join(indir, run2d, 'fieldlist-'+run2d+'.fits')

    if field is not None and mjd is not None:
        if epoch is True:
            field_dir = ptt.join(indir, run2d, field, 'epoch')
        else:
            field_dir = ptt.join(indir, run2d, field)
    elif (custom is not None) and (mjd is not None):
        field_dir = ptt.join(indir, run2d, custom)
    
    
    if logfile is None:
        if outroot is not None:
            logfile = ptt.join(outroot+'.log')
        elif (field is not None) and (mjd is not None):
            logfile = ptt.join(field_dir,'spAll-'+field+'-'+mjd+'.log')
        elif (custom is not None) and (mjd is not None):
            logfile = ptt.join(field_dir,'spAll-'+custom+'-'+mjd+'.log')
        elif field is None and mjd is None:
            if epoch is True:
                logfile = ptt.join(indir, run2d,'spAll_epoch-'+run2d+'.log')
            elif custom is None:
                logfile = ptt.join(indir, run2d,'spAll-'+run2d+'.log')
            else:
                logfile = ptt.join(indir, f'spAll_{custom}-{run2d}.log')
        else:
            if epoch is True:
                logfile = ptt.join(indir, 'spAll_epoch.log')
            elif custom is None:
                logfile = ptt.join(indir, 'spAll.log')
            else:
                logfile = ptt.join(indir, f'spAll_{custom}.log')
        if dev:
            logfile = logfile.replace('spAll','spAll_dev')
        
    splog.open(logfile=logfile, backup=False)
    splog.log('Log file '+logfile+' opened '+ time.ctime())
    
    cflist = False
    if allsky is False:
        flist = None
        if ptt.exists(fieldlist_file):
            splog.log(f'Reading fieldlist file: {fieldlist_file}')
            try:
                flist = Table(fits.getdata(fieldlist_file))
            except:
                time.sleep(90)
                try:
                    flist = Table(fits.getdata(fieldlist_file))
                except:
                    pass
        else:
            time.sleep(30)
            if ptt.exists(fieldlist_file):
                splog.log(f'Reading fieldlist file: {fieldlist_file}')
                try:
                    flist = Table(fits.getdata(fieldlist_file))
                except:
                    time.sleep(90)
                    try:
                        flist = Table(fits.getdata(fieldlist_file))
                    except:
                        pass
        if (flist is None) and (dev is False):
            cflist = True
            if field is not None and mjd is not None:
                fmlog = f'fieldlist-{field}-{mjd}.log'
            else:
                fmlog = None
            splog.log('ERROR: '+ptt.basename(fieldlist_file)+' is required (Running Now)')
            flist = fieldlist(create=True, topdir=indir, run2d=[run2d], run1d=[getenv('RUN1D')],
                             outdir=None, legacy=legacy, custom=custom, basehtml=None, epoch=epoch,
                             logfile=fmlog, noplot=True)
    else:
        flist = build_custom_fieldlist(indir, custom, run2d, run1d)

    full_flist = flist    
    #if run2d is not None:
    #    flist = flist[np.where(flist['RUN2D'] == run2d)[0]]
    if programs is not None:
        ftemp = None
        for prog in programs:
            if ftemp is None:
                ftemp = flist[np.where(flist['PROGRAM'] == prog)[0]]
            else:
                ftemp = vstack([ftemp, flist[np.where(flist['PROGRAM'] == prog)[0]]])
        flist = ftemp

    if allsky is False:
        if (field is not None) and (mjd is not None):
            flist = flist[np.logical_and((flist['FIELD'].data.astype(int) == int(field)), (flist['MJD'].data.astype(int) == int(mjd)))]
    else:
        if mjd is not None:
            flist = flist[(flist['MJD'] == int(mjd))]
    spAll = None
    spline = None
    spAll_fmjds = None
    spline_fmjds = None

    if (field is not None) and (mjd is not None):
        idx  = np.where((flist['FIELD'] == field) & (flist['MJD'] == int(mjd)))[0]
        if len(idx) == 0:
            flist.add_row({'FIELD': field, 'MJD':mjd, 'STATUS2D': 'unknown'})

    spallfile, spalllitefile, splinefile, spAlldatfile = build_fname(indir, run2d, outroot=outroot,
                                                                     field=field, mjd=mjd, dev=dev,
                                                                     epoch=epoch, allsky=allsky,
                                                                     custom=custom)
    if not clobber:
        if ptt.exists(spallfile):
            splog.log(f'Reading Existing spAll file: {spallfile}')
            spAll = Table(fits.getdata(spallfile))
            try:
                spAll_fmjds = spAll['FIELD','MJD']
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD'])
            except:
                spAll_fmjds = Table(names = ['FIELD','MJD'])
        elif ptt.exists(spallfile.replace('.gz','')):
            splog.log(f"Reading Existing spAll file: {spallfile.replace('.gz','')}")
            spAll = Table(fits.getdata(spallfile.replace('.gz','')))
            try:
                spAll_fmjds = spAll['FIELD','MJD']
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD'])
            except:
                spAll_fmjds = Table(names = ['FIELD','MJD'])
        else:
            spAll_fmjds = Table(names = ['FIELD','MJD'])
        if ptt.exists(splinefile):
            splog.log(f'Reading Existing spLine file: {splinefile}')
            spline = Table(fits.getdata(splinefile))
            try:
                spline_fmjds = spline['FIELD','MJD']
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD'])
            except:
                spAll_fmjds = Table(names = ['FIELD','MJD'])
        elif ptt.exists(splinefile.replace('.gz','')):
            splog.log(f"Reading Existing spLine file: {splinefile.replace('.gz','')}")
            spline = Table(fits.getdata(splinefile.replace('.gz','')))
            try:
                spline_fmjds = spline['FIELD','MJD']
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD'])
            except:
                spAll_fmjds = Table(names = ['FIELD','MJD'])
        else:
            spline_fmjds = Table(names = ['FIELD','MJD'])
        
        if remerge_fmjd is not None:
            idx  = np.where((spAll['FIELD'] == int(remerge_fmjd.split('-')[0])) & (spAll['MJD'] == int(remerge_fmjd.split('-')[1])))[0]
            idxl = np.where((spline['FIELD'] == int(remerge_fmjd.split('-')[0])) & (spline['MJD'] == int(remerge_fmjd.split('-')[1])))[0]
            if len(idx) > 0:
                spAll.remove_rows(idx)
                spAll_fmjds = spAll['FIELD','MJD']
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD'])
            if len(idxl) > 0:
                spline.remove_rows(idxl)
                spline_fmjds = spline['FIELD','MJD']
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD'])
    flist.sort(['MJD','FIELD'])
    j = 0
    for i, row in enumerate(flist):
        if spAll_fmjds is not None:
                
            idx = spAll_fmjds[(spAll_fmjds['FIELD'] == int(row['FIELD'])) & (spAll_fmjds['MJD'] == int(row['MJD']))]
            idxl = spline_fmjds[(spline_fmjds['FIELD'] == int(row['FIELD'])) & (spline_fmjds['MJD'] == int(row['MJD']))]
            if len(idx)*len(idxl) > 0:
                splog.log(f"Skipping (Complete) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
                continue
        if (field is not None) and (mjd is not None):
            if (row['STATUS2D'].lower().strip() != 'done') or (row['STATUSCOMBINE'].lower().strip() != 'done') or (row['STATUS1D'].lower().strip() != 'done'):
                splog.log(f"Checking incomplete status ({row['STATUS2D'].strip()} RUN2D) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
                fmlog = f'fieldlist-{field}-{mjd}.log'
                row = fieldlist(create=True, topdir=indir, run2d=[run2d], run1d=[getenv('RUN1D')],
                                  outdir=None, legacy=legacy, custom=custom, basehtml=None, epoch=epoch,
                                  logfile=fmlog, field=field, mjd=mjd, noplot=True, fmsplog=splog)
                flist[i] = row
        if row['STATUS2D'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUS2D'].strip()} RUN2D) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
            continue
        elif row['STATUSCOMBINE'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUSCOMBINE'].strip()} Combine) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
            continue
        elif row['STATUS1D'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUS1D'].strip()} RUN1D) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
            continue
        
        splog.log(f"Reading/Building Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
        onefield = oneField(row, row['FIELD'], row['MJD'],
                            skip_line=skip_line, include_bad=include_bad,
                            legacy=legacy, skip_specprimary=skip_specprimary, dev=dev,
                            XCSAO=XCSAO, indir=indir, clobber=clobber, epoch = epoch,
                            merge_only=merge_only, custom = custom, allsky = allsky)
        if not merge_only:
            write_spAll(onefield['spall'], onefield['spline'], None, indir, run2d, datamodel,
                        line_datamodel, outroot=None, field=row['FIELD'], mjd = row['MJD'],
                        verbose=verbose, dev=dev, clobber=clobber, epoch=epoch, silent=True,
                        custom = custom, allsky = allsky, SDSSC2BV = SDSSC2BV)
        if onefield['spall'] is None:
            continue
        if not skip_line:
            if onefield['spline'] is None:
                continue
                
        if spAll is None:
            spAll = onefield['spall']
        else:
            try:
                spAll = vstack([spAll,onefield['spall']])
            except:
                for col in spAll.colnames:
                    if len(spAll[col].shape) <=1: continue
                    if spAll[col].shape[1] > onefield['spall'][col].shape[1]:
                        coldat = onefield['spall'][col].data
                        pad = spAll[col].shape[1] - onefield['spall'][col].shape[1]
                        onefield['spall'][col] = np.pad(coldat, [(0,0),(pad,0)], mode = 'constant', constant_values= 0)
                    elif spAll[col].shape[1] < onefield['spall'][col].shape[1]:
                        coldat = spAll[col].data
                        pad = onefield['spall'][col].shape[1] - spAll[col].shape[1]
                        spAll[col] = np.pad(coldat, [(0,0),(pad,0)], mode = 'constant', constant_values= 0)
                spAll = vstack([spAll,onefield['spall']])
            
        if onefield['spline'] is not None:
            if spline is None:
                spline = onefield['spline']
            else:
                spline = vstack([spline,onefield['spline']])
        if  merge_only:
            if limit is not None:
                j = j+1
                del onefield
                gc.collect()
                splog.log(f'{j}:{limit} ({time.ctime()})')
                if j >= limit:
                    break
    if custom is not None and mjd is not None: field = custom
    if not(field is not None and mjd is not None):
        if spline is not None:
            if len(spline) == 0:
                spline = None
        if allsky is True:
            spall_raw = spAll.copy()
            spAll.sort('MJD')
            spAll = unique(spAll, keys='CATALOGID', keep='last')
            spAll.sort('CATALOGID')
            if spline is not None:
                dropped = setdiff(spall_raw, spAll, keys=['TARGET_INDEX','MJD','FIELD'])
                dropped = dropped['TARGET_INDEX','MJD','FIELD','CATALOGID']
                if len(dropped) > 0:
                    dropped = unique(dropped)
                    for row in dropped:

                        idx = np.where((spline['TARGET_INDEX'].data == row['TARGET_INDEX']) &
                                        (spline['MJD'].data == row['MJD']) &
                                        (spline['FIELD'].data == row['FIELD']))[0]
                                    
                        spline.remove_rows(idx)
        if spAll is None:
            splog.info('No valid spAll entries')
            splog.info('EXITING!!')
            exit()
        if not skip_specprimary:
            spAll = specPrimary(spAll)

        if lite is True:
            spAll_lite = spAll.copy()
            for i in range(len(spAll)):
                with warnings.catch_warnings():
                    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                    try:
                        spAll_lite[i]['ASSIGNED']            = str(min(np.array(spAll[i]['ASSIGNED'].split()).astype(int)))
                    except:
                        spAll_lite[i]['ASSIGNED']            = '0'
                    try:
                        spAll_lite[i]['ON_TARGET']           = str(min(np.array(spAll[i]['ON_TARGET'].split()).astype(int)))
                    except:
                        spAll_lite[i]['ON_TARGET']           = '0'
                    try:
                        spAll_lite[i]['VALID']               = str(min(np.array(spAll[i]['VALID'].split()).astype(int)))
                    except:
                        spAll_lite[i]['VALID']               = '0'
                    try:
                        spAll_lite[i]['DECOLLIDED']          = str(min(np.array(spAll[i]['DECOLLIDED'].split()).astype(int)))
                    except:
                        spAll_lite[i]['DECOLLIDED']          = '0'
                    try:
                        spAll_lite[i]['MOON_DIST']           = str(np.nanmean(np.array(spAll[i]['MOON_DIST'].split()).astype(float)))
                    except:
                        spAll_lite[i]['MOON_DIST']           = 'nan'
                    try:
                        spAll_lite[i]['MOON_PHASE']          = str(np.nanmean(np.array(spAll[i]['MOON_PHASE'].split()).astype(float)))
                    except:
                        spAll_lite[i]['MOON_PHASE']          = 'nan'
                    try:
                        spAll_lite[i]['CARTON_TO_TARGET_PK'] = str(int(np.array(spAll[i]['CARTON_TO_TARGET_PK'].split())[0]))
                    except:
                        spAll_lite[i]['CARTON_TO_TARGET_PK'] = '0'
            
            for col in ['ASSIGNED','ON_TARGET','VALID','DECOLLIDED']:
                spAll_lite[col].fill_value = False
                spAll_lite[col] = spAll_lite[col].astype(bool)

            for col in ['CARTON_TO_TARGET_PK']:
                spAll_lite[col] = spAll_lite[col].astype(int)

            for col in ['MOON_DIST','MOON_PHASE']:
                spAll_lite[col] = spAll_lite[col].astype(float)
        else:
            spAll_lite = None

        write_spAll(spAll, spline, spAll_lite, indir, run2d, datamodel, line_datamodel,
                    epoch=epoch, dev=dev, outroot=outroot, field=field, mjd=mjd,
                    verbose=verbose, clobber=True, custom = custom, allsky = allsky,
                    SDSSC2BV = SDSSC2BV)

    if spAll is None:
        splog.info('No valid spAll entries')
        splog.info('EXITING!!')
        exit()

    if allsky is False:
        #flist = Table(fits.getdata(fieldlist_file))
        plot_sky(ptt.dirname(fieldlist_file), full_flist, fieldlist_file)

    if field is not None and mjd is not None:
        splog.log(f'Successful completion of fieldmerge for {field}-{mjd} at '+ time.ctime())
    elif custom is not None:
        splog.log(f'Successful completion of fieldmerge for {custom}-{mjd} at '+ time.ctime())
    else:
        splog.log('Successful completion of fieldmerge at '+ time.ctime())
    splog.close()
    return

def build_fname(indir, run2d, outroot=None, field=None, mjd=None, dev=False,
                epoch=False, custom=None, allsky=False):
    if outroot is not None:
        spallfile     = ptt.join(outroot+'.fits.gz')
        spalllitefile = ptt.join(outroot+'-lite'+'.fits.gz')
        splinefile    = ptt.join(outroot+'Line'+'.fits.gz')
        spAlldatfile  = ptt.join(outroot+'.dat.gz')
    else:
        if custom is not None:
            if allsky is True:
                field = custom
            elif field is not None:
                field = field_to_string(field)
        elif field is not None:
            field = field_to_string(field)
        if field is not None and mjd is not None:
            field = field_to_string(field)
            mjd = str(mjd)
            if epoch is True:
                specfull_dir  = ptt.join(indir, run2d, 'epoch', 'spectra','full', field, mjd)
            else:
                specfull_dir  = ptt.join(indir, run2d, 'spectra','full', field, mjd)
            spallfile     = ptt.join(specfull_dir, 'spAll-'+field+'-'+mjd+'.fits.gz')
            spalllitefile = ptt.join(specfull_dir, 'spAll-lite-'+field+'-'+mjd+'.fits.gz')
            splinefile    = ptt.join(specfull_dir, 'spAllLine-'+field+'-'+mjd+'.fits.gz')
            spAlldatfile  = ptt.join(specfull_dir, 'spAll-'+field+'-'+mjd+'.dat.gz')
        elif run2d is not None:
            if epoch is True:
                spAll_dir  = ptt.join(indir, run2d, 'epoch')
            else:
                spAll_dir  = ptt.join(indir, run2d)
                
            if custom is None:
                spallfile     = ptt.join(spAll_dir, 'spAll-'+run2d+'.fits.gz')
                spalllitefile = ptt.join(spAll_dir, 'spAll-lite-'+run2d+'.fits.gz')
                splinefile    = ptt.join(spAll_dir, 'spAllLine-'+run2d+'.fits.gz')
                spAlldatfile  = ptt.join(spAll_dir, 'spAll-'+run2d+'.dat.gz')
            else:
                spallfile     = ptt.join(spAll_dir, 'spAll-'+run2d+'-'+custom+'.fits.gz')
                spalllitefile = ptt.join(spAll_dir, 'spAll-lite-'+run2d+'-'+custom+'.fits.gz')
                splinefile    = ptt.join(spAll_dir, 'spAllLine-'+run2d+'-'+custom+'.fits.gz')
                spAlldatfile  = ptt.join(spAll_dir, 'spAll-'+run2d+'-'+custom+'.dat.gz')
        else:
            if epoch is True:
                spAll_dir  = ptt.join(indir, 'epoch')
            else:
                spAll_dir  = ptt.join(indir)
            spallfile     = ptt.join(spAll_dir, 'spAll.fits.gz')
            spalllitefile = ptt.join(spAll_dir, 'spAll-lite.fits.gz')
            splinefile    = ptt.join(spAll_dir, 'spAllLine.fits.gz')
            spAlldatfile  = ptt.join(spAll_dir, 'spAll.dat.gz')
    if dev:
        spallfile = spallfile.replace('spAll','spAll_dev')
        spalllitefile = spalllitefile.replace('spAll','spAll_dev')
        splinefile = splinefile.replace('spAllLine','spAllLine_dev')
        spAlldatfile = spAlldatfile.replace('spAll','spAll_dev')
    
    return(spallfile, spalllitefile, splinefile, spAlldatfile)

def write_spAll(spAll, spline, spAll_lite, indir, run2d, datamodel, line_datamodel,
                epoch=False, dev=False, outroot=None, field=None, mjd = None,
                verbose=False, clobber=False, silent=False, custom = None,
                allsky = False, SDSSC2BV = ''):
 
    spallfile, spalllitefile, splinefile, spAlldatfile = build_fname(indir, run2d, outroot=outroot,
                                                                     field=field, mjd=mjd, dev=dev,
                                                                     epoch=epoch, custom=custom, allsky=allsky)

    drop_cols = None
    date = time.ctime()
    if spline is not None:
        spline = merge_dm(table=spline, ext = 'spZline', name = 'SPLINE', dm = line_datamodel,
                          splog=splog, drop_cols=drop_cols, verbose=verbose)
    if spAll_lite is not None:
        #spall_lite = spall_lite[np.where(spall_lite['CATALOGID'] != -999)[0]]
        spAll_lite = merge_dm(table=spAll_lite, ext = 'SPALL_lite', name = 'SPALL', dm = datamodel,
                              splog=splog,drop_cols=drop_cols, verbose=verbose)
    exists = ptt.exists(spallfile) if not clobber else False
    if spAll is None: return
    if not exists:
        hdul = merge_dm(ext = 'Primary', hdr = {'RUN2D':run2d,
                                                'Date':time.ctime(),
                                                'SDSSC2BV': SDSSC2BV},
                        dm = datamodel, splog=splog, verbose=verbose)
        spAll = merge_dm(table=spAll, ext = 'SPALL', name = 'SPALL', dm = datamodel, splog=splog,drop_cols=drop_cols, verbose=verbose)
        makedirs(ptt.dirname(spallfile), exist_ok=True)
        splog.log('Writing '+spallfile)
        fits.HDUList([hdul,spAll]).writeto(spallfile, overwrite=True, checksum=True)
        del hdul, spAll
        gc.collect()
    elif not silent:
        splog.log('Skipping '+spallfile+' (exists)')
    if spAll_lite is not None:
        exists = ptt.exists(spalllitefile) if not clobber else False
        if not exists:
            makedirs(ptt.dirname(spalllitefile), exist_ok=True)
            hdul_lite = merge_dm(ext = 'Primary', hdr = {'RUN2D':run2d,
                                                         'Date':time.ctime(),
                                                         'SDSSC2BV': SDSSC2BV},
                                 dm = datamodel, splog=splog, verbose=verbose)
            fits.HDUList([hdul_lite, spAll_lite]).writeto(spalllitefile, overwrite=True, checksum=True)
            splog.log('Writing '+spalllitefile)
            del hdul_lite, spAll_lite
            gc.collect()
        elif not silent:
            splog.log('Skipping '+spalllitefile+' (exists)')
    if spline is not None:
        exists = ptt.exists(splinefile) if not clobber else False
        if not exists:
            makedirs(ptt.dirname(splinefile), exist_ok=True)
            hdul_line = merge_dm(ext = 'Primary', hdr = {'RUN2D':run2d,'Date':time.ctime()},
                                 dm = line_datamodel, splog=splog, verbose=verbose)
            fits.HDUList([hdul_line,spline]).writeto(splinefile, overwrite=True, checksum=True)
            splog.log('Writing '+splinefile)
            del hdul_line, spline
            gc.collect()
        elif not silent:
            splog.log('Skipping '+splinefile+' (exists)')

    #spAll.write(spAlldatfile, format='ascii.fixed_width_two_line')
    
    return
    
if __name__ == '__main__' :
    """ 
    build spAll    
    """
    parser = argparse.ArgumentParser(
            prog=ptt.basename(sys.argv[0]),
            description='Build BOSS spAll Summary File')

    parser.add_argument('--run2d', type=str, help='Optional override value for the enviro variable $RUN2D', default=getenv('RUN2D'))
    parser.add_argument('--indir', type=str, help='Optional override value for the environment variable $BOSS_SPECTRO_REDUX', default = getenv('BOSS_SPECTRO_REDUX'))
    

    parser.add_argument('--skip_line',        action='store_true', help='skip the generation of spAllLine.fits')
    parser.add_argument('--include_bad',      action='store_true', help='include bad fields')
    parser.add_argument('--legacy',           action='store_true', help='Include columns used by SDSS-IV and depreciated in SDSS-V')
    parser.add_argument('--skip_specprimary', action='store_true', help='Skip creation of specprimary and associated columns')
    parser.add_argument('--lite',             action='store_true', help='Produce lite version of spAll file')
    parser.add_argument('--XCSAO',            action='store_true', help='Include XCSAO columns')
    parser.add_argument('--field', '-f',      type=str,            help='Run for a single Field', default=None)
    parser.add_argument('--mjd',   '-m',      type=str,            help='Run for a single MJD', default=None)
    parser.add_argument('--clobber',          action='store_true', help='Clobber all spAll-field-mjd files')
    parser.add_argument('--verbose',          action='store_true', help='Log columns not saved')
    parser.add_argument('--logfile',          type=str,            help='Manually set logfile')
    parser.add_argument('--epoch',            action='store_true', help='Produce spAll for epoch coadds')
    parser.add_argument('--dev',              action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--programs',         nargs='*',           help='List of programs to include')
    parser.add_argument('--datamodel',        type=str,            help='Supply a spAll datamodel file (defaults to $IDLSPEC2D/datamodel/spall_dm.par')
    parser.add_argument('--line_datamodel',   type=str,            help='Supply a spline datamodel file (defaults to $IDLSPEC2D/datamodel/spzline_dm.par')
    parser.add_argument('--outroot',          type=str,            help='Path and root of filename for output (defaults to $BOSS_SPECTRO_REDUX/$RUN2D/{field}/{mjd}/spAll)')
    parser.add_argument('--remerge_fmjd', '-r',type=str,           help='Field-MJD to replace in spAll')
    parser.add_argument('--merge_only', '-o', action='store_true', help='Skip Building new spAll-Field-MJD files and just merge existing')
    parser.add_argument('--allsky',           action='store_true', help='Build spAll for Allsky Custom Coadd')
    parser.add_argument('--custom',           type=str,            help='Name of Custom Coadd')
    parser.add_argument('--run1d',            type=str,            help='Optional override value for the enviro variable $RUN1D (only for custom allsky coadds)', default=getenv('RUN1D'))
    parser.add_argument('--limit',            type=int,            help='Limit number of Field-MJD spAll files to read before save', default = None)

    args = parser.parse_args()
    
    if args.merge_only is True:
        args.clobber = False
    fieldmerge(**vars(args))
