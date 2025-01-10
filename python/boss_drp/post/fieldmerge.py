#!/usr/bin/env python3
from boss_drp.summary import Summary_names, summary_names, fieldlist_name
from boss_drp.post.fieldlist import fieldlist
from boss_drp.field import field_to_string, Field, fieldgroup
from boss_drp.utils import (merge_dm, get_lastline)
from boss_drp.utils import match as wwhere
from boss_drp.post import plot_sky_targets, plot_sky_locations
from boss_drp.utils import specobjid, retry
from boss_drp.utils.splog import splog

from sdss_semaphore.targeting import TargetingFlags
import traceback
import argparse
import sys
import os.path as ptt
from os import getenv, makedirs, remove, rename
import numpy as np
from astropy.io import fits
from astropy.table import (Table, Column, vstack, join,
                          unique, setdiff, hstack,MaskedColumn)
from healpy import ang2pix
from pydl.pydlutils.spheregroup import spheregroup
import time
from datetime import timedelta, datetime
import warnings
from glob import glob
import gc
import shutil

run2d_warn = True


def read_zans(field_class):
    zansfile = ptt.join(field_class.spec1d_dir(),
                    'spZbest-' + field_class.field + '-' + field_class.mjd + '.fits')
    if ptt.exists(zansfile):
        splog.log('Reading Zbest file: '+ptt.basename(zansfile))
        try:
            zans = Table.read(zansfile)
            for key in zans.colnames:
                zans.rename_column(key,key.upper())
            zans.meta = {}
        except:
            splog.log('Error reading '+zansfile)
            zans = None
    else:
        splog.log(zansfile +' missing')
        zans = None
    return(zans)

def read_zline(field_class):
    zlinefile = ptt.join(field_class.spec1d_dir(),
                    'spZline-' + field_class.field + '-' + field_class.mjd + '.fits')
    if ptt.exists(zlinefile):
        splog.log('Reading Zline file: '+ptt.basename(zlinefile))
        try:
            zline = Table.read(zlinefile)
            for key in zline.colnames:
                zline.rename_column(key,key.upper())
            zline.add_column(field_class.obs, name = 'OBS')
            zline.meta = {}
        except:
            splog.log('Error reading '+zlinefile)
            zline = None
    else:
        splog.log(zlinefile +' missing')
        zline = None
    return(zline)

def read_xcsao(field_class):#sp1d_dir, field, mjd, spAll):
    XCSAOfile = ptt.join(field_class.spec1d_dir(),
                'spXCSAO-' + field_class.field + '-' + str(field_class.mjd) + '.fits')
    if ptt.exists(XCSAOfile):
        splog.log('Reading XCSAO file: '+ptt.basename(XCSAOfile))
        try:
            XCSAO_tab = Table.read(XCSAOfile,1)
            for key in XCSAO_tab.colnames:
                XCSAO_tab.rename_column(key,key.upper())
            XCSAO_tab.meta = {}
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

def read_fibermap(field_class, allsky=False):
    spfieldfile = ptt.join(field_class.dir(),
            'spField-'+field_class.field+'-'+field_class.mjd+'.fits')
    if allsky is True:
        spfieldfile = spfieldfile.replace('spField-','spFullsky-')
    if ptt.exists(spfieldfile):
        splog.log('Reading Fibermap: '+ptt.basename(spfieldfile))
        try:
            fibermap = Table.read(spfieldfile, 5)
            fibermap.meta = {}
            hdr = fits.getheader(spfieldfile,0)
            fibermap['SPEC_FILE'] = np.char.add(np.char.add('spec-'+field_class.field+'-'+field_class.mjd+'-',
                                        np.char.strip(np.asarray(fibermap['CATALOGID']).astype(str))),'.fits')
            renames = {'XFOCAL_LIST':'XFOCAL','YFOCAL_LIST':'YFOCAL',
                        'CARTON_TO_TARGET_PK_LIST':'CARTON_TO_TARGET_PK',
                        'ASSIGNED_LIST':'ASSIGNED', 'ON_TARGET_LIST':'ON_TARGET',
                        'VALID_LIST':'VALID','DECOLLIDED_LIST':'DECOLLIDED',
                        'TOO_LIST':'TOO','ICATALOGID':'CATALOGID',
                        'GAIA_ID_DR2': 'GAIA_ID', 'MJDLIST': 'MJD_LIST',
                        'PROGRAM':'PROGRAMNAME'}
        
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

            if 'CATALOGID_V0' not in fibermap.colnames:
                fibermap.add_column(-999, name = 'CATALOGID_V0')
                fibermap.add_column(-999, name = 'CATALOGID_V0P5')
                fibermap.add_column(-999, name = 'CATALOGID_V1')
                catid = fibermap['CATALOGID']
                cv0   = fibermap['CATALOGID_V0']
                cv0p5 = fibermap['CATALOGID_V0P5']
                cv1   = fibermap['CATALOGID_V1']
                iv0   = np.where(wwhere(fibermap['CATVERSION'].data.astype(str),'0.0*'))[0]
                iv0p1 = np.where(wwhere(fibermap['CATVERSION'].data.astype(str),'0.1*'))[0]
                iv0p5 = np.where(wwhere(fibermap['CATVERSION'].data.astype(str),'0.5*'))[0]
                iv1   = np.where(wwhere(fibermap['CATVERSION'].data.astype(str),'1.*'))[0]
                cv0[iv0]     = catid[iv0]
                cv0[iv0p1]   = catid[iv0p1]
                cv0p5[iv0p5] = catid[iv0p5]
                cv1[iv1]     = cv1[iv1]

            HEALPIX_PATH = np.char.add(np.char.add('$MWM_HEALPIX/', np.asarray(fibermap['HEALPIXGRP']).astype(str)),'/')
            HEALPIX_PATH = np.char.add(np.char.add(HEALPIX_PATH, np.asarray(fibermap['HEALPIX']).astype(str)),'/boss/')
            HEALPIX_PATH = np.char.add(HEALPIX_PATH, fibermap['SPEC_FILE'])
            idx = np.where(fibermap['CATALOGID'] == 0)[0]
            HEALPIX_PATH[idx] = ''
            fibermap.add_column(HEALPIX_PATH, name = 'HEALPIX_PATH')
        except:
            splog.log('Error reading '+spfieldfile)
            fibermap = Table()
    else:
        splog.log(spfieldfile +' missing')
        fibermap = Table()
    return(fibermap)

def read_spline(splinefile, skip_line=False, clobber=False, merge_only=False):
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
            spline = Table.read(splinefile)
            if len(spline) == 0:
                splog.debug('Empty spline file: '+ptt.basename(splinefile))
                remove(splinefile)
                spline = None
            spline.meta = {}
        except:
            splog.log('Failure opening '+ splinefile)
            spline = None
    elif merge_only:
        spline = None
        splog.info(f'No Existing spline file ({ptt.basename(splinefile)})')
    else:
        spline = None
    return(spline)

def read_spall(spAllfile, clobber=False, merge_only=False):
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
            splog.info('Reading spAll file: '+ptt.basename(spAllfile))
            spAll = Table.read(spAllfile)
            if len(spAll) == 0:
                splog.debug('Empty spAll file: '+ptt.basename(spAllfile))
                remove(spAllfile)
                spAll = None
            spAll.meta = {}
        except:
            splog.info('Failure opening '+ spAllfile)
            spAll = None
    elif merge_only:
        spAll = None
        splog.info(f'No Existing spAll file ({ptt.basename(spAllfile)})')
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
    
    field_class = Field(indir,row['RUN2D'], field, custom_name=custom,
                        epoch =epoch, run1d=row['RUN1D'],
                        mjd = mjd, obs=row['OBSERVATORY'])
    field_dir = field_class.dir()
    sp1d_dir = field_class.spec1d_dir(row['RUN1D'])
    dev = False
    fnames = Summary_names()
    fnames.set(indir, row['RUN2D'], field=field, mjd=mjd, dev=dev,
               epoch=epoch, custom=custom, allsky=allsky)
    
    spAll  = read_spall(fnames.spAllfile, clobber=clobber, merge_only=merge_only)
    spline = read_spline(fnames.splinefile, skip_line = skip_line, clobber=clobber, merge_only=merge_only)

    if merge_only:
        return({'spall':spAll, 'spline': spline}, fnames)


    ############################
    #  Build SpAll and spAllline
    ############################
    zans = None
    if spAll is None:
        zans = read_zans(field_class)
        spAll = read_fibermap(field_class,allsky=allsky)#field_dir, field, mjd, epoch=epoch, allsky=allsky)
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
                XCSAO_tab = read_xcsao(field_class)
                if XCSAO_tab is not None:
                    spAll = join(spAll,XCSAO_tab,keys='TARGET_INDEX', join_type='left')

            if zans is None:
                return ({'spall':None, 'spline': None}, fnames)
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
                spAll = build_specobjid(spAll, epoch = epoch, custom = custom)
                
    if spline is None and skip_line is False:
        spline = read_zline(field_class)#sp1d_dir, field, mjd, row['OBSERVATORY'])

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
    return({'spall':spAll, 'spline': spline}, fnames)



def build_specobjid(spAll,custom=None, epoch=False):
    splog.info('Building SPECOBJIDs')
    if epoch:
        coadd = 'epoch'
    elif custom is not None:
        coadd = custom
    else:
        coadd = 'daily'
    spAll.add_column(Column('',name='SPECOBJID', dtype=object))
    spAll['SPECOBJID'] = specobjid.encode(spAll['SDSS_ID'].data,
                                          spAll['FIELD'].data.astype(str),
                                          spAll['MJD'].data, coadd,
                                          spAll['RUN2D'].data.astype(str),
                                          fiberid=spAll['TARGET_INDEX'].data,
                                          allnew = False)

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

def specPrimary_sdssid(spAll, update = False):
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
            
    specprim = spAll['SPECPRIMARY']
    specboss = spAll['SPECBOSS']
    nspecobs = spAll['NSPECOBS']
    if isinstance(spAll['SDSS_ID'], MaskedColumn):
        sdssids = spAll['SDSS_ID'].filled(-999).data
    else:
        sdssids = spAll['SDSS_ID'].data
    sidx = np.where((sdssids < 0))[0]
    specprim[sidx] = -999
    specboss[sidx] = -999
    nspecobs[sidx] = -999

    if update:
        splog.log(f'Keeping current and updating only')
        idx_new = np.where(specprim == -1)[0]
        sdssids = np.unique(sdssids[idx_new])
    else:
        spAll['SPECPRIMARY'] = 0
        spAll['SPECBOSS'] = 0
        spAll['NSPECOBS'] = 0
        spAll['BOSS_SPECOBJ_ID'] = 0

    sdssids = np.unique(sdssids)
    sdssids = sdssids[sdssids >= 0]
    counter_step = round(len(sdssids)/100)
    for i, id in enumerate(sdssids):
        if (i % counter_step) == 0:
            if i + counter_step < len(sdssids):
                splog.info(f'Assigning SpecPrimaries: {i+1} - {i+100000} (of {len(sdssids)})')
            else:
                splog.info(f'Assigning SpecPrimaries: {i+1} - {len(sdssids)} (of {len(sdssids)})')
        idx = np.where(spAll['SDSS_ID'].data == id)[0]
        nspecobs[idx] = len(idx)
        if update:
            specprim[idx] = 0
            specboss[idx] = 0
        primary = np.argmax(score[idx])
        specprim[idx[primary]] = 1
        specboss[idx[primary]] = 1
    splog.log('Time to assign primaries = '+str(timedelta(seconds = time.time() - t2)))
    return(spAll)


def build_custom_fieldlist(indir, custom, run2d, run1d):
    flist = Table()
    spfields = []
    for cc in [custom,custom+'_lco',custom+'_apo']:
        fc = Field(indir, run2d, custom_name = cc, custom = True)
        spfields.extend(glob(ptt.join(fc.dir(), f'spFullsky-{cc}-*.fits')))
    for spfield in spfields:
        hdr = fits.getheader(spfield,0)
        
        cc = ptt.basename(spfield).split('-')[1]
        if cc.split('_')[-1] in ['lco','apo']:
            obs = cc.split('_')[-1].upper()
        else:
            obs = ''
        fc = Field(indir, run2d, custom_name = cc, custom = True)
        spdiagcomblog = ptt.join(fc.dir(),f"spDiagcomb-{cc}-{str(hdr['RUNMJD'])}.log")
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
        spDiag1dlog = ptt.join(fc.spec1d_dir(run1d), f"spDiag1d-{cc}-{str(hdr['MJD'])}.log")
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
            
        if obs.lower() == 'lco':
            sn2 = [sn2_g2,sn2_i2]
        elif obs.lower() == 'apo':
            sn2 = [sn2_g1,sn2_i1]
        else:
            sn2 = [sn2_g1,sn2_i1,sn2_g2,sn2_i2]
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')

            flist_row = Table({'RUN2D':[run2d], 'RUN1D':[run1d],
                               'PROGRAMNAME':[custom], 'FIELD':[0], 'MJD': [hdr['MJD']],
                               'STATUS2D':['Done'], 'STATUSCOMBINE':[STATUSCOMBINE],
                               'STATUS1D':[STATUS1D],'FIELDQUALITY':['good'],
                               'FIELDSN2':[np.nanmin( np.array(sn2))],
                               'OBSERVATORY':[obs]})

        flist =  vstack([flist, flist_row])
    flist.pprint_all()
    return(flist)


def fieldmerge(run2d=getenv('RUN2D'), indir= getenv('BOSS_SPECTRO_REDUX'),
               skip_line=False, include_bad=False, legacy=False, skip_specprimary=False,
               update_specprimary=True, lite=False, XCSAO=False, field=None,
               mjd=None, programs=None, clobber=False, dev=False,
               datamodel=None, line_datamodel=None, verbose=False, epoch =False, outroot=None,
               logfile=None, remerge_fmjd=None, remerge_mjd=None, merge_only=False, limit=None,
               custom=None, allsky=False, run1d=None, bkup=False, mjdstart=None):
               
    try:
        SDSSC2BV = TargetingFlags.meta['SDSSC2BV']
    except:
        SDSSC2BV = '1'
    if field is not None:
        field = field_to_string(field)
    if mjd is not None:
        mjd = str(mjd)
        
    if datamodel is not None:
        summary_names.datamodel = datamodel
    if line_datamodel is not None:
        summary_names.line_datamodel = line_datamodel

    fieldlist_name.build(indir, run2d, epoch=epoch, custom_name=custom)
    fieldlist_file = fieldlist_name.name
    if field is not None and mjd is not None:
        fc = Field(indir, run2d, field, custom_name = custom, epoch = epoch)
        field_dir = fc.dir()
    
    
    if logfile is None:
        if outroot is not None:
            logfile = ptt.join(outroot+'.log')
        elif (field is not None) and (mjd is not None):
            logfile = ptt.join(field_dir,'spAll-'+field+'-'+mjd+'.log')
        elif (custom is not None) and (mjd is not None):
            logfile = ptt.join(field_dir,'spAll-'+custom+'-'+mjd+'.log')
        elif field is None and mjd is None:
            if epoch is True:
                logfile = ptt.join(indir, run2d,'summary','epoch','spAll_epoch-'+run2d+'.log')
            elif custom is None:
                logfile = ptt.join(indir, run2d,'summary','daily','spAll-'+run2d+'.log')
            else:
                logfile = ptt.join(indir, run2d,'summary',fieldgroup(custom, custom=True),f'spAll_{custom}-{run2d}.log')
        else:
            if epoch is True:
                logfile = ptt.join(indir, 'spAll_epoch.log')
            elif custom is None:
                logfile = ptt.join(indir, 'spAll.log')
            else:
                logfile = ptt.join(indir, f'spAll_{custom}.log')
        if dev:
            logfile = logfile.replace('spAll','spAll_dev')
    if len(ptt.dirname(logfile)) > 0:
        makedirs(ptt.dirname(logfile),exist_ok=True)
    splog.open(logfile=logfile, backup=False)
    splog.log('Log file '+logfile+' opened '+ time.ctime())
    
    cflist = False
    if allsky is False:
        flist = None
        if ptt.exists(fieldlist_file):
            splog.log(f'Reading fieldlist file: {fieldlist_file}')
            try:
                flist = Table.read(fieldlist_file)
            except:
                time.sleep(90)
                try:
                    flist = Table(fieldlist_file)
                except:
                    pass
        else:
            time.sleep(30)
            if ptt.exists(fieldlist_file):
                splog.log(f'Reading fieldlist file: {fieldlist_file}')
                try:
                    flist = Table.read(fieldlist_file)
                except:
                    time.sleep(90)
                    try:
                        flist = Table.read(fieldlist_file)
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
                             logfile=fmlog, noplot=True, return_tab = True)
    else:
        flist = build_custom_fieldlist(indir, custom, run2d, run1d)

    flist = flist['RUN2D','RUN1D','PROGRAMNAME','FIELD','MJD','STATUS2D',
                  'STATUSCOMBINE','STATUS1D','FIELDQUALITY','FIELDSN2','OBSERVATORY']
    if programs is not None:
        ftemp = None
        for prog in programs:
            if ftemp is None:
                ftemp = flist[np.where(flist['PROGRAMNAME'] == prog)[0]]
            else:
                ftemp = vstack([ftemp, flist[np.where(flist['PROGRAMNAME'] == prog)[0]]])
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

    
    if (field is not None) and (mjd is not None) and (custom is None):
        idx  = np.where((flist['FIELD'] == field) & (flist['MJD'] == int(mjd)))[0]
        if len(idx) == 0:
            flist.add_row({'FIELD': field, 'MJD':mjd, 'STATUS2D': 'unknown'})

    summary_names.set(indir, run2d, outroot=outroot, field=field, mjd=mjd,
                      dev = dev, epoch=epoch, allsky=allsky, custom=custom)
    if not clobber:
        spAll = None
        spline = None
        if ptt.exists(summary_names.spAllfile):
            splog.log(f'Reading Existing spAll file: {summary_names.spAllfile}')
            spAll = Table.read(summary_names.spAllfile)
        elif ptt.exists(summary_names.spAllfile.replace('.gz','')):
            splog.log(f"Reading Existing spAll file: {summary_names.spAllfile.replace('.gz','')}")
            try:
                spAll = Table.read(summary_names.spAllfile.replace('.gz',''))
            except:
                time.sleep(60)
                spAll = Table.read(summary_names.spAllfile.replace('.gz',''))
        else:
            spAll_fmjds = Table(names = ['FIELD','MJD','OBS'])
        if spAll is not None:
            try:
                spAll_fmjds = spAll['FIELD','MJD','OBS']
                
                if isinstance(spAll_fmjds['FIELD'], MaskedColumn):
                    spAll_fmjds['FIELD'] = spAll_fmjds['FIELD'].filled(0)
                if isinstance(spAll_fmjds['OBS'], MaskedColumn):
                    spAll_fmjds['OBS'] = spAll_fmjds['OBS'].filled('???')
                if isinstance(spAll_fmjds['MJD'], MaskedColumn):
                    spAll_fmjds['MJD'] = spAll_fmjds['MJD'].filled(0)
                    
                try:
                    spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD','OBS'])
                except Exception as e:
                    splog.warning(f'{type(e).__name__}: {e}')
                    splog.info('Rerunning unique')
                    spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD'])
            except Exception as e:
                splog.warning(f'{type(e).__name__}: {e}')
                spAll_fmjds = Table(names = ['FIELD','MJD','OBS'])
        try:
            spAll.meta = {}
            spAll_fmjds.meta = {}
        except:
            pass
        if ptt.exists(summary_names.splinefile):
            splog.log(f'Reading Existing spLine file: {summary_names.splinefile}')
            try:
                spline = Table.read(summary_names.splinefile)
            except:
                time.sleep(60)
                spline = Table.read(summary_names.splinefile)
        elif ptt.exists(summary_names.splinefile.replace('.gz','')):
            splog.log(f"Reading Existing spLine file: {summary_names.splinefile.replace('.gz','')}")
            try:
                spline = Table.read(splinefile.replace('.gz',''))
            except:
                time.sleep(60)
                spline = Table.read(splinefile.replace('.gz',''))
        else:
            spline_fmjds = Table(names = ['FIELD','MJD','OBS'])
        if spline is not None:
            try:
                try:
                    spline_fmjds = spline['FIELD','MJD','OBS']
                except:
                    spline_fmjds = spline['FIELD','MJD']
                    spline_fmjds.add_column('???', name='OBS')
                    
                if isinstance(spline_fmjds['FIELD'], MaskedColumn):
                    spline_fmjds['FIELD'] = spline_fmjds['FIELD'].filled(0)
                if isinstance(spline_fmjds['OBS'], MaskedColumn):
                    spline_fmjds['OBS'] = spline_fmjds['OBS'].filled('???')
                if isinstance(spline_fmjds['MJD'], MaskedColumn):
                    spline_fmjds['MJD'] = spline_fmjds['MJD'].filled(0)
                    
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD','OBS'])
                
            except Exception as e:
                splog.warning(f'{type(e).__name__}: {e}')
                spline_fmjds = Table(names = ['FIELD','MJD','OBS'])
        try:
            spline.meta = {}
            spline_fmjds.meta = {}
        except:
            pass
        if remerge_fmjd is not None:
            idx  = np.where((spAll_fmjds['FIELD'] == int(remerge_fmjd.split('-')[0])) &
                            (spAll_fmjds['MJD'] == int(remerge_fmjd.split('-')[1])))[0]
            idxl = np.where((spline_fmjds['FIELD'] == int(remerge_fmjd.split('-')[0])) &
                            (spline_fmjds['MJD'] == int(remerge_fmjd.split('-')[1])))[0]
            if len(idx) > 0:
                spAll_fmjds.remove_rows(idx)
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD','OBS'])
            if len(idxl) > 0:
                spline_fmjds.remove_rows(idxl)
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD','OBS'])
        if remerge_mjd is not None:
            idx  = np.where((spAll_fmjds['MJD']  == int(remerge_mjd)))[0]
            idxl = np.where((spline_fmjds['MJD'] == int(remerge_mjd)))[0]
            if len(idx) > 0:
                spAll_fmjds.remove_rows(idx)
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD','OBS'])
            if len(idxl) > 0:
                spline_fmjds.remove_rows(idxl)
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD','OBS'])
    flist.sort(['MJD','FIELD'])
    if mjdstart is not None:
        splog.info(f"Only Checking Field-MJDs with MJD >= {mjdstart}")
        flist = flist[flist['MJD'] >= mjdstart]
    j = 0
    for i, row in enumerate(flist):
        tfield_str = f"Field:{row['FIELD']}  MJD:{row['MJD']} OBS:{row['OBSERVATORY']} ({i+1}/{len(flist)})"
        if spAll_fmjds is not None:
                
            idx = spAll_fmjds[(spAll_fmjds['FIELD'] == int(row['FIELD'])) &
                              (spAll_fmjds['MJD'] == int(row['MJD'])) &
                              (spAll_fmjds['OBS'] == row['OBSERVATORY'])]
                    
            idxl = spline_fmjds[(spline_fmjds['FIELD'] == int(row['FIELD'])) &
                                (spline_fmjds['MJD'] == int(row['MJD'])) &
                                ((spline_fmjds['OBS'] == row['OBSERVATORY']) |
                                 (spline_fmjds['OBS'] == '???'))]

            if len(idx)*len(idxl) > 0:
                splog.log(f"Skipping (Complete) {tfield_str}")
                continue
        if (field is not None) and (mjd is not None) and (allsky is False):
            if (row['STATUS2D'].lower().strip() != 'done') or (row['STATUSCOMBINE'].lower().strip() != 'done') or (row['STATUS1D'].lower().strip() != 'done'):
                splog.log(f"Checking incomplete status ({row['STATUS2D'].strip()} RUN2D) {tfield_str}")
                fmlog = f'fieldlist-{field}-{mjd}.log'
                row = retry(fieldlist, retries=3, delay = 5, logger=splog.log,
                            create=True, topdir=indir, run2d=[run2d], run1d=[getenv('RUN1D')],
                            outdir=None, legacy=legacy, custom=custom, basehtml=None, epoch=epoch,
                            logfile=fmlog, field=field, mjd=mjd, noplot=True)
                    
                try:
                    flist[i] = row
                except:
                    try:
                        flist[i] = row[flist.colnames]
                    except:

                        print(flist.colnames)
                        print('-----------------')
                        print(row)
                        print('------------------')
                        print(row[flist.colnames])

                        print('------------------')
                        print(len(row))
                        print(len(flist.columns))
                        raise
        if row['STATUS2D'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUS2D'].strip()} RUN2D) {tfield_str}")
            continue
        elif row['STATUSCOMBINE'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUSCOMBINE'].strip()} Combine) {tfield_str}")
            continue
        elif row['STATUS1D'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUS1D'].strip()} RUN1D) {tfield_str}")
            continue
        
        splog.log(f"Reading/Building {tfield_str}")
        if dev:
            dev1 = True if merge_only else False
        else:
            dev1 = False
        field_clobber = clobber if not merge_only else False
        if custom is not None:
            rfield = '_'.join([row['PROGRAMNAME'],row['OBSERVATORY'].lower()])
        else:
            rfield = row['FIELD']
        onefield, fnames = oneField(row, rfield, row['MJD'],
                                    skip_line=skip_line, include_bad=include_bad,
                                    legacy=legacy, skip_specprimary=skip_specprimary,
                                    dev=dev1, XCSAO=XCSAO, indir=indir,
                                    clobber=field_clobber, epoch = epoch,
                                    merge_only=merge_only, custom = custom, allsky = allsky)
        if onefield['spall'] is None:
            continue
        if not merge_only:
            fnames.datamodel = summary_names.datamodel
            fnames.line_datamodel = summary_names.line_datamodel
            write_spAll(onefield['spall'].copy(), onefield['spline'].copy(), None,
                        run2d, fnames, verbose=verbose, clobber=clobber, silent=True,
                        SDSSC2BV = SDSSC2BV)

        onefield['spall']['SPECPRIMARY'] = -1
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
                        if col.upper() == 'SDSS5_TARGET_FLAGS':
                            onefield['spall'][col] = np.pad(coldat, [(0,0),(0,pad)],
                                                            mode = 'constant',
                                                            constant_values= 0)
                        else:
                            onefield['spall'][col] = np.pad(coldat, [(0,0),(pad,0)],
                                                            mode = 'constant',
                                                            constant_values= 0)
                        coldat = None
                        pad = None
                    elif spAll[col].shape[1] < onefield['spall'][col].shape[1]:
                        coldat = spAll[col].data
                        pad = onefield['spall'][col].shape[1] - spAll[col].shape[1]
                        if col.upper() == 'SDSS5_TARGET_FLAGS':
                            spAll[col] = np.pad(coldat, [(0,0),(0,pad)],
                                                mode = 'constant',
                                                constant_values= 0)
                        else:
                            spAll[col] = np.pad(coldat, [(0,0),(pad,0)],
                                                mode = 'constant',
                                                constant_values= 0)
                        coldat = None
                        pad = None
                spAll = vstack([spAll,onefield['spall']])
            
        if onefield['spline'] is not None:
            if spline is None:
                spline = onefield['spline']
            else:
                spline = vstack([spline,onefield['spline']])
        if (merge_only) and (limit is not None):
            j = j+1
            del onefield
            gc.collect()
            splog.log(f'{j}:{limit} ({time.ctime()})')
            if j >= limit:
                break
        else:
            del onefield
            gc.collect()
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
                dropped = setdiff(spall_raw, spAll, keys=['TARGET_INDEX','MJD','OBS'])
                dropped = dropped['TARGET_INDEX','MJD','OBS','CATALOGID']
                if len(dropped) > 0:
                    dropped = unique(dropped)
                    for row in dropped:

                        idx = np.where((spline['TARGET_INDEX'].data == row['TARGET_INDEX']) &
                                        (spline['MJD'].data == row['MJD']))[0]
                                    
                        spline.remove_rows(idx)
        if spAll is None:
            splog.info('No valid spAll entries')
            splog.info('EXITING!!')
            exit()
        if not skip_specprimary:
            #spAll = specPrimary(spAll)
            spAll = specPrimary_sdssid(spAll, update=update_specprimary)

        if lite is True:
            splog.info('Creating spAll-lite Table')
            mr = len(spAll)



            for col in ['ASSIGNED','ON_TARGET','VALID','DECOLLIDED', 'TOO','CARTON_TO_TARGET_PK']:
                if col not in spAll.columns:
                    splog.info(f'{col} missing from spAll')
                    spAll[col] = '-999'
            for col in ['MOON_DIST','MOON_PHASE','DELTA_RA_LIST','DELTA_DEC_LIST']:
                if col not in spAll.columns:
                    splog.info(f'{col} missing from spAll')
                    spall[col] = 'nan'

            spAll_lite = spAll['ASSIGNED','ON_TARGET','VALID','DECOLLIDED', 'TOO',
                               'MOON_DIST','MOON_PHASE','CARTON_TO_TARGET_PK',
                               'DELTA_RA_LIST','DELTA_DEC_LIST'].copy()
            errors = {}
            for i in range(mr):
                if (i % 100000) == 0:
                    if i + 100000 < mr:
                        splog.info(f'Re-Formatting arrays in spAll-lite rows: {i+1} - {i+100000} (of {mr})')
                    else:
                        splog.info(f'Re-Formatting arrays in spAll-lite rows: {i+1} - {mr} (of {mr})')
                        
                for col in ['ASSIGNED','ON_TARGET','VALID','DECOLLIDED','TOO','CARTON_TO_TARGET_PK']:
                    try:
                        if col in ['CARTON_TO_TARGET_PK']:
                            spAll_lite[col][i] = str(int(np.asarray(spAll[col][i].split())[0]))
                        else:
                            try:
                                spAll_lite[col][i] = str(min(np.asarray(spAll[col][i].split()).astype(int)))
                            except:
                                spAll_lite[col][i] = str(min(np.asarray(spAll[col][i].split()).astype(float).astype(int)))
                    except Exception as e:
                        if f'{col}: {type(e).__name__}: {e}' not in errors.keys():
                            errors[f'{col}: {type(e).__name__}: {e}'] = 1
                            splog.warning(f'{col}: {type(e).__name__}: {e}')
                        else:
                            errors[f'{col}: {type(e).__name__}: {e}'] += 1
                            
                        tb = traceback.extract_tb(e.__traceback__)
                        filename, line, func, text = tb[-1]
                        print(f'{filename}:{func}:{line}:{col}: {type(e).__name__}: {e}', file=sys.stderr)
                        spAll_lite[col][i] = '0'
                with warnings.catch_warnings():
                    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                    for col in ['MOON_DIST','MOON_PHASE','DELTA_RA_LIST','DELTA_DEC_LIST']:
                        try:
                            temp = np.nanmean(np.asarray(spAll[col][i].split()).astype(float))
                            column_dtype = spAll_lite[col].dtype
                            if column_dtype.kind == 'U':  # Unicode string
                                prec = column_dtype.itemsize // 4  # 4 bytes per character
                            elif column_dtype.kind == 'S':  # Byte string
                                prec = column_dtype.itemsize
                            else: # unknown format, just concert and accept the warning if it truncates
                                spAll_lite[col][i] = str(temp)
                                continue
                            temps = "{:.{prec}g}".format(temp, prec=prec-1)
                            if '.' not in temps and 'e' not in temps:
                                temps += ".0" # if there is not trailing zero
                            if '.' in temps and 'e' not in temps:
                                temps = temps[:prec] # if the string is to long (due to floating point errors) then trim it
                                temps = temps.rstrip('0').rstrip('.') # remove excess trailing zeros
                            spAll_lite[col][i] = temps

                        except Exception as e:
                            filename, line, func, text = tb[-1]
                            print(f'{filename}:{func}:{line}:{col}: {type(e).__name__}: {e}', file=sys.stderr)
                            if f'{col}: {type(e).__name__}: {e}' not in errors.keys():
                                errors[f'{col}: {type(e).__name__}: {e}'] = 1
                                splog.warning(f'{col}: {type(e).__name__}: {e}')
                            else:
                                errors[f'{col}: {type(e).__name__}: {e}'] += 1
                            spAll_lite[col][i] = 'nan'

            if len(errors) > 0:
                splog.warning('----------------------\n spAll->spAll-liteConverstion Error Summary\n----------------------')
                for er in errors:
                    splog.warning(f'{er}: {errors[er]}rows')
                splog.warning('----------------------')
            for col in ['ASSIGNED','ON_TARGET','VALID','DECOLLIDED','TOO','CARTON_TO_TARGET_PK']:
                spAll_lite[col].fill_value = -999
                try:
                    spAll_lite[col] = spAll_lite[col].astype(int)
                except Exception as e:
                    splog.warning(f'{type(e).__name__}: {e}')
                    column_data = np.array(spAll_lite[col])
                    column_data = np.array([x.decode('utf-8') if x is not None else None for x in column_data])
                    column_data = np.asarray([x == 'True' if x is not None else None for x in column_data]).astype(int)
                    spAll_lite.remove_column(col)
                    spAll_lite[col] = MaskedColumn(column_data, mask=[x is None for x in column_data])
                    column_data = None
                    del column_data

            spAll_lite.rename_column('DELTA_RA_LIST', 'DELTA_RA')
            spAll_lite.rename_column('DELTA_DEC_LIST', 'DELTA_DEC')

            for col in ['CARTON_TO_TARGET_PK']:
                spAll_lite[col] = spAll_lite[col].astype(int)

            for col in ['MOON_DIST','MOON_PHASE','DELTA_RA','DELTA_DEC']:
                try:
                    spAll_lite[col] = spAll_lite[col].astype(float)
                except Exception as e:
                    splog.warning(f'{type(e).__name__}: {e} - {col}')
                    spAll_lite[col].set_fill_value(np.NaN)
                    spAll_lite[col] = spAll_lite[col].astype(float)

        else:
            spAll_lite = None

        write_spAll(spAll, spline, spAll_lite, run2d,
                                summary_names, verbose=verbose, clobber=True,
                                SDSSC2BV = SDSSC2BV, bkup=bkup)

    if spAll is None:
        splog.info('No valid spAll entries')
        splog.info('EXITING!!')
        exit()

    if allsky is False:
        summary_names.set(indir, run2d, outroot=outroot, dev=dev, epoch=epoch,
                               custom=custom, allsky=allsky)



        spAll = None
        spline = None
        spAll_lite = None
        if not(field is not None and mjd is not None):
            plot_sky_locations()
            plot_sky_targets(nobs=True)

    if field is not None and mjd is not None:
        splog.log(f'Successful completion of fieldmerge for {field}-{mjd} at '+ time.ctime())
    elif custom is not None:
        splog.log(f'Successful completion of fieldmerge for {custom}-{mjd} at '+ time.ctime())
    else:
        splog.log('Successful completion of fieldmerge at '+ time.ctime())
    splog.close()
    return



def write_spAll(spAll, spline, spAll_lite, run2d, fnames,
                verbose=False, clobber=False, silent=False,
                SDSSC2BV = '', bkup = False):

    drop_cols = None
    date = time.ctime()
    exists = ptt.exists(fnames.spAllfile) if not clobber else False
    if spAll is None: return
    if bkup:
        if ptt.exists(fnames.spAllfile):
            try:
                bkup_str = datetime.strptime(fits.getheader(fnames.spAllfile,0)['DATE'],'%c').isoformat()
            except:
                bkup_str = (datetime.now() - timedelta(days=1)).isoformat()
            fnames.bk.set(flag = bkup_str)
            if not ptt.exists(fnames.bk.spAllfile):
                shutil.copy2(fnames.spAllfile,fnames.bk.spAllfile)
                shutil.copy2(fnames.spAlllitefile,fnames.bk.spAlllitefile)
                shutil.copy2(fnames.splinefile,fnames.bk.splinefile)

    
    if not exists:
        hdul = merge_dm(ext = 'Primary', hdr = {'RUN2D':run2d,
                                                'Date':time.ctime(),
                                                'SDSSC2BV': SDSSC2BV},
                        dm = fnames.datamodel, verbose=verbose)
        splog.info('Formatting spAll table')
        spAll = merge_dm(table=spAll, ext = 'SPALL', name = 'SPALL', dm = fnames.datamodel,
                         drop_cols=drop_cols, verbose=verbose)
        makedirs(ptt.dirname(fnames.spAllfile), exist_ok=True)
        splog.log('Writing '+fnames.spAllfile)
        fits.HDUList([hdul,spAll]).writeto(fnames.temp.spAllfile,
                                           overwrite=True, checksum=True)
        hdul = None
        if spAll_lite is not None:
            spAll = Table(spAll.data)
            for col in spAll_lite.colnames:
                try:
                    spAll.remove_column(col)
                except Exception as e:
                    if col in ['DELTA_RA', 'DELTA_DEC']:
                        continue
                    splog.warning(f'{col} {type(e).__name__}: {e}')
                    pass
            spAll = hstack([spAll, spAll_lite], join_type = 'exact')
        else:
            spAll = None
        gc.collect()
    elif not silent:
        splog.log('Skipping '+fnames.spAllfile+' (exists)')
    if spAll_lite is not None:
        exists = ptt.exists(fnames.spAlllitefile) if not clobber else False
        if not exists:
            makedirs(ptt.dirname(fnames.spAlllitefile), exist_ok=True)
            hdul_lite = merge_dm(ext = 'Primary', hdr = {'RUN2D':run2d,
                                                         'Date':time.ctime(),
                                                         'SDSSC2BV': SDSSC2BV},
                                 dm = fnames.datamodel,verbose=verbose)
            splog.info('Formatting spAll-lite table')
            spAll = merge_dm(table=spAll, ext = 'SPALL_lite', name = 'SPALL', dm = fnames.datamodel,
                                  drop_cols=drop_cols, verbose=verbose)
            splog.log('Writing '+fnames.spAlllitefile)
            fits.HDUList([hdul_lite, spAll]).writeto(fnames.temp.spAlllitefile,
                                                          overwrite=True, checksum=True)
            del hdul_lite, spAll_lite
            gc.collect()
        elif not silent:
            splog.log('Skipping '+fnames.spAlllitefile+' (exists)')
    if spline is not None:
        exists = ptt.exists(fnames.splinefile) if not clobber else False
        if not exists:
            makedirs(ptt.dirname(fnames.splinefile), exist_ok=True)
            hdul_line = merge_dm(ext = 'Primary', hdr = {'RUN2D':run2d,'Date':time.ctime()},
                                 dm = fnames.line_datamodel, verbose=verbose)
            splog.info('Formatting spAllLine table')
            spline = merge_dm(table=spline, ext = 'spZline', name = 'SPLINE', dm = fnames.line_datamodel,
                              drop_cols=drop_cols, verbose=verbose)
            splog.log('Writing '+fnames.splinefile)
            fits.HDUList([hdul_line,spline]).writeto(fnames.temp.splinefile,
                                                     overwrite=True, checksum=True)
            del hdul_line, spline
            gc.collect()
        elif not silent:
            splog.log('Skipping '+fnames.splinefile+' (exists)')

    try:
        rename(fnames.temp.spAllfile,fnames.spAllfile)
    except:
        pass
    try:
        rename(fnames.temp.spAlllitefile,fnames.spAlllitefile)
    except:
        pass
    try:
        rename(fnames.temp.splinefile,fnames.splinefile)
    except:
        pass

    #spAll.write(spAlldatfile, format='ascii.fixed_width_two_line')
    
    return
    
