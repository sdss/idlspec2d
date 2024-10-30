#!/usr/bin/env python3
from boss_drp import idlspec2d_dir
from boss_drp.post.fieldlist import fieldlist
from boss_drp.field import field_to_string, field_spec_dir, fieldgroup
from boss_drp.field import field_dir as create_field_dir
from boss_drp.utils import (merge_dm, Splog, get_lastline)
from boss_drp.utils import match as wwhere
from boss_drp.post import plot_sky_targets, plot_sky_locations
from boss_drp.utils import specobjid, retry

from sdss_semaphore.targeting import TargetingFlags

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
splog = Splog()


def read_zans(sp1d_dir, field, mjd):
    zansfile = ptt.join(sp1d_dir,'spZbest-' + field + '-' + mjd + '.fits')
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

def read_zline(sp1d_dir, field, mjd):
    zlinefile = ptt.join(sp1d_dir,'spZline-' + field + '-' + mjd + '.fits')
    if ptt.exists(zlinefile):
        splog.log('Reading Zline file: '+ptt.basename(zlinefile))
        try:
            zline = Table.read(zlinefile)
            for key in zline.colnames:
                zline.rename_column(key,key.upper())
            zline.meta = {}
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
            fibermap = Table.read(spfieldfile, 5)
            fibermap.meta = {}
            hdr = fits.getheader(spfieldfile,0)
            fibermap['SPEC_FILE'] = np.char.add(np.char.add('spec-'+field+'-'+mjd+'-',
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
            spline = Table.read(splinefile)
            if len(spline) == 0:
                remove(splinefile)
                spline = None
            spline.meta = {}
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
            spAll = Table.read(spAllfile)
            if len(spAll) == 0:
                remove(spAllfile)
                spAll = None
            spAll.meta = {}
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
            field_dir =  create_field_dir(ptt.join(indir,row['RUN2D']), field)
            sp1d_dir =  ptt.join(field_dir, 'epoch', row['RUN1D'])
        else:
            field_dir =  create_field_dir(ptt.join(indir,row['RUN2D']), field)
            sp1d_dir =  ptt.join(field_dir, row['RUN1D'])
    else:
        if epoch is True:
            field_dir =  create_field_dir(ptt.join(indir,row['RUN2D']), field, custom=True)
            sp1d_dir =  ptt.join(field_dir, 'epoch', row['RUN1D'])
        else:
            field_dir =  create_field_dir(ptt.join(indir,row['RUN2D']), field, custom=True)
            sp1d_dir =  ptt.join(field_dir, row['RUN1D'])
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
                spAll = build_specobjid(spAll, epoch = epoch, custom = custom)
                
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
    if update:
        splog.log(f'Keeping current and updating only')
        specprim = spAll['SPECPRIMARY']
        idx_new = np.where(specprim == -1)[0]
        try:
            sdssids = spAll['SDSS_ID'].filled(-999).data
        except:
            sdssids = spAll['SDSS_ID'].data
        sdssids = np.unique(sdssids[idx_new])
        sdssids = sdssids[sdssids >= 0]
    else:
        try:
            sdssids = np.unique(spAll['SDSS_ID'].filled(-999).data)
        except:
            sdssids = np.unique(spAll['SDSS_ID'].data)
        sdssids = sdssids[sdssids >= 0]
        spAll['SPECPRIMARY'] = 0
        spAll['SPECBOSS'] = 0
        spAll['NSPECOBS'] = 0
    for id in sdssids:
        idx = np.where(spAll['SDSS_ID'].data == id)[0]
        spAll[idx]['NSPECOBS'] = len(idx)
        if update:
            spAll[idx]['SPECPRIMARY'] = 0
            spAll[idx]['SPECBOSS'] = 0
        primary = np.argmax(score[idx])
        spAll[idx[primary]]['SPECPRIMARY'] = 1
        spAll[idx[primary]]['SPECBOSS'] = 1
    splog.log('Time to assign primaries = '+str(timedelta(seconds = time.time() - t2)))
    return(spAll)


def build_custom_fieldlist(indir, custom, run2d, run1d):
    flist = Table()
    topdir2d = ptt.join(indir, run2d)
    spfields = []
    for cc in [custom,custom+'_lco',custom+'_apo']:
        spfields.extend(glob(ptt.join(create_field_dir(topdir2d,cc, custom=True),
                                 f'spFullsky-{cc}-*.fits')))
    for spfield in spfields:
        hdr = fits.getheader(spfield,0)
        
        cc = ptt.basename(spfield).split('-')[1]
        if cc.split('_')[-1] in ['lco','apo']:
            obs = cc.split('_')[-1].upper()
        else:
            obs = ''
        spdiagcomblog = ptt.join(create_field_dir(topdir2d,cc, custom=True),
                                 f"spDiagcomb-{cc}-{str(hdr['RUNMJD'])}.log")
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
        spDiag1dlog = ptt.join(create_field_dir(topdir2d,cc, custom=True),
                               run1d, f"spDiag1d-{cc}-{str(hdr['MJD'])}.log")
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
                           'PROGRAMNAME':[custom], 'FIELD':[0], 'MJD': [hdr['MJD']],
                           'STATUS2D':['Done'], 'STATUSCOMBINE':[STATUSCOMBINE],
                           'STATUS1D':[STATUS1D],'FIELDQUALITY':['good'],
                           'FIELDSN2':[np.nanmin( np.array([sn2_g1,sn2_i1,sn2_g2,sn2_i2]))],
                           'OBSERVATORY':[obs]})
    
        flist =  vstack([flist, flist_row])
    print(flist)
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
        
    if datamodel is None:
        datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spall_dm.par')
    if line_datamodel is None:
        line_datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spzline_dm.par')
 
    dflags = [indir, run2d,'summary']
    if custom is not None:
        dflags.append(custom)
        fflag = '-'+custom
    elif epoch is True:
        dflags.append('epoch')
        fflag = '-epoch'
    else:
        dflags.append('daily')
        fflag =''
    fieldlist_file = ptt.join(*dflags, 'fieldlist-'+run2d+fflag+'.fits')

    if (custom is not None) and (mjd is not None):
#        if field is not None:
#            custom = field
#            field = None
        field_dir = create_field_dir(ptt.join(indir, run2d), field, custom=True)
    elif field is not None and mjd is not None:
        field_dir = create_field_dir(ptt.join(indir, run2d),field)
        if epoch is True:
            field_dir = ptt.join(field_dir, 'epoch')
    
    
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

    spallfile, spalllitefile, splinefile, spAlldatfile = build_fname(indir, run2d, outroot=outroot,
                                                                     field=field, mjd=mjd, dev=dev,
                                                                     epoch=epoch, allsky=allsky,
                                                                     custom=custom)
    if not clobber:
        if ptt.exists(spallfile):
            splog.log(f'Reading Existing spAll file: {spallfile}')
            spAll = Table.read(spallfile)
            try:
                spAll_fmjds = spAll['FIELD','MJD']
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD'])
            except:
                spAll_fmjds = Table(names = ['FIELD','MJD'])
        elif ptt.exists(spallfile.replace('.gz','')):
            splog.log(f"Reading Existing spAll file: {spallfile.replace('.gz','')}")
            try:
                spAll = Table.read(spallfile.replace('.gz',''))
            except:
                time.sleep(60)
                spAll = Table.read(spallfile.replace('.gz',''))
            try:
                spAll_fmjds = spAll['FIELD','MJD']
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD'])
            except:
                spAll_fmjds = Table(names = ['FIELD','MJD'])
        else:
            spAll_fmjds = Table(names = ['FIELD','MJD'])
        try:
            spAll.meta = {}
            spAll_fmjds.meta = {}
        except:
            pass
        if ptt.exists(splinefile):
            splog.log(f'Reading Existing spLine file: {splinefile}')
            try:
                spline = Table.read(splinefile)
            except:
                time.sleep(60)
                spline = Table.read(splinefile)
            try:
                spline_fmjds = spline['FIELD','MJD']
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD'])
            except:
                spline_fmjds = Table(names = ['FIELD','MJD'])
        elif ptt.exists(splinefile.replace('.gz','')):
            splog.log(f"Reading Existing spLine file: {splinefile.replace('.gz','')}")
            try:
                spline = Table.read(splinefile.replace('.gz',''))
            except:
                time.sleep(60)
                spline = Table.read(splinefile.replace('.gz',''))
            try:
                spline_fmjds = spline['FIELD','MJD']
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD'])
            except:
                spline_fmjds = Table(names = ['FIELD','MJD'])
        else:
            spline_fmjds = Table(names = ['FIELD','MJD'])
        try:
            spline.meta = {}
            spline_fmjds.meta = {}
        except:
            pass
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
        if remerge_mjd is not None:
            idx  = np.where((spAll['MJD']  == int(remerge_mjd)))[0]
            idxl = np.where((spline['MJD'] == int(remerge_mjd)))[0]
            if len(idx) > 0:
                spAll.remove_rows(idx)
                spAll_fmjds = spAll['FIELD','MJD']
                spAll_fmjds = unique(spAll_fmjds,keys=['FIELD','MJD'])
            if len(idxl) > 0:
                spline.remove_rows(idxl)
                spline_fmjds = spline['FIELD','MJD']
                spline_fmjds = unique(spline_fmjds,keys=['FIELD','MJD'])
    flist.sort(['MJD','FIELD'])
    if mjdstart is not None:
        splog.info(f"Only Checking Field-MJDs with MJD >= {mjdstart}")
        flist = flist[flist['MJD'] >= mjdstart]
    j = 0
    for i, row in enumerate(flist):
        if spAll_fmjds is not None:
                
            idx = spAll_fmjds[(spAll_fmjds['FIELD'] == int(row['FIELD'])) & (spAll_fmjds['MJD'] == int(row['MJD']))]
            idxl = spline_fmjds[(spline_fmjds['FIELD'] == int(row['FIELD'])) & (spline_fmjds['MJD'] == int(row['MJD']))]
            if len(idx)*len(idxl) > 0:
                splog.log(f"Skipping (Complete) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
                continue
        if (field is not None) and (mjd is not None) and (allsky is False):
            if (row['STATUS2D'].lower().strip() != 'done') or (row['STATUSCOMBINE'].lower().strip() != 'done') or (row['STATUS1D'].lower().strip() != 'done'):
                splog.log(f"Checking incomplete status ({row['STATUS2D'].strip()} RUN2D) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
                fmlog = f'fieldlist-{field}-{mjd}.log'
                row = retry(fieldlist, retries=3, delay = 5, logger=splog.log,
                            create=True, topdir=indir, run2d=[run2d], run1d=[getenv('RUN1D')],
                            outdir=None, legacy=legacy, custom=custom, basehtml=None, epoch=epoch,
                            logfile=fmlog, field=field, mjd=mjd, noplot=True, fmsplog=splog)
                    
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
            splog.log(f"Skipping ({row['STATUS2D'].strip()} RUN2D) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
            continue
        elif row['STATUSCOMBINE'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUSCOMBINE'].strip()} Combine) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
            continue
        elif row['STATUS1D'].lower().strip() != 'done':
            splog.log(f"Skipping ({row['STATUS1D'].strip()} RUN1D) Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
            continue
        
        splog.log(f"Reading/Building Field:{row['FIELD']}  MJD:{row['MJD']} ({i+1}/{len(flist)})")
        
        if dev:
            dev1 = True if merge_only else False
        else:
            dev1 = False
        field_clobber = clobber if not merge_only else False
        if custom is not None:
            rfield = '_'.join([row['PROGRAMNAME'],row['OBSERVATORY'].lower()])
        else:
            rfield = row['FIELD']
        onefield = oneField(row, rfield, row['MJD'],
                            skip_line=skip_line, include_bad=include_bad,
                            legacy=legacy, skip_specprimary=skip_specprimary, dev=dev1,
                            XCSAO=XCSAO, indir=indir, clobber=field_clobber, epoch = epoch,
                            merge_only=merge_only, custom = custom, allsky = allsky)
        if onefield['spall'] is None:
            continue
        if not merge_only:
            write_spAll(onefield['spall'].copy(), onefield['spline'].copy(), None,
                        indir, run2d, datamodel,
                        line_datamodel, outroot=None, field=rfield, mjd = row['MJD'],
                        verbose=verbose, dev=dev, clobber=clobber, epoch=epoch, silent=True,
                        custom = custom, allsky = allsky, SDSSC2BV = SDSSC2BV)

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
                                        (spline['MJD'].data == row['MJD']) &
                                        (spline['OBS'].data == row['OBS']))[0]
                                    
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



            for col in ['ASSIGNED','ON_TARGET','VALID','DECOLLIDED', 'TOO',
                        'CARTON_TO_TARGET_PK']:
                if col not in spAll.columns:
                    spAll[col] = '0'
            for col in ['MOON_DIST','MOON_PHASE']:
                if col not in spAll.columns:
                    spall[col] = 'nan'

            spAll_lite = spAll['ASSIGNED','ON_TARGET','VALID','DECOLLIDED', 'TOO',
                               'MOON_DIST','MOON_PHASE','CARTON_TO_TARGET_PK'].copy()
            for i in range(mr):
                if (i % 100000) == 0:
                    if i + 100000 < mr:
                        splog.info(f'Re-Formatting arrays in spAll-lite rows: {i+1} - {i+100000} (of {mr})')
                    else:
                        splog.info(f'Re-Formatting arrays in spAll-lite rows: {i+1} - {mr} (of {mr})')
                with warnings.catch_warnings():
                    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                    for col in ['ASSIGNED','ON_TARGET','VALID','DECOLLIDED','TOO']:
                        try:
                            spAll_lite[i][col]            = str(min(np.array(spAll[i][col].split()).astype(int)))
                        except:
                            spAll_lite[i][col]            = '0'
                    for col in ['MOON_DIST','MOON_PHASE']:
                        try:
                            spAll_lite[i][col]           = str(np.nanmean(np.array(spAll[i][col].split()).astype(float)))
                        except:
                            spAll_lite[i][col]           = 'nan'
                    for col in ['CARTON_TO_TARGET_PK']:
                        try:
                            spAll_lite[i][col] = str(int(np.array(spAll[i][col].split())[0]))
                        except:
                            spAll_lite[i][col] = '0'
            for col in ['ASSIGNED','ON_TARGET','VALID','DECOLLIDED','TOO']:
                spAll_lite[col].fill_value = False
                try:
                    spAll_lite[col] = spAll_lite[col].astype(bool)
                except:
                    column_data = np.array(spAll_lite[col])
                    column_data = np.array([x.decode('utf-8') if x is not None else None for x in column_data])
                    column_data = np.array([x == 'True' if x is not None else None for x in column_data])
                    spAll_lite.remove_column(col)
                    spAll_lite[col] = MaskedColumn(column_data, mask=[x is None for x in column_data])
                    column_data = None
                    del column_data


            for col in ['CARTON_TO_TARGET_PK']:
                spAll_lite[col] = spAll_lite[col].astype(int)

            for col in ['MOON_DIST','MOON_PHASE']:
                spAll_lite[col] = spAll_lite[col].astype(float)


        else:
            spAll_lite = None

        spAll_file = write_spAll(spAll, spline, spAll_lite, indir, run2d, datamodel, line_datamodel,
                    epoch=epoch, dev=dev, outroot=outroot, field=field, mjd=mjd,
                    verbose=verbose, clobber=True, custom = custom, allsky = allsky,
                    SDSSC2BV = SDSSC2BV, bkup=bkup)

    if spAll is None:
        splog.info('No valid spAll entries')
        splog.info('EXITING!!')
        exit()

    if allsky is False:
        spallfile, spalllitefile, splinefile, spAlldatfile = build_fname(indir, run2d, outroot=outroot, dev=dev,
                                                                        epoch=epoch, custom=custom, allsky=allsky)



        spAll = None
        spline = None
        spAll_lite = None
        if not(field is not None and mjd is not None):
            outdir = ptt.dirname(fieldlist_file)
            fieldlist_file = ptt.basename(fieldlist_file)
            plot_sky_locations(outdir, fieldlist_file, splog)
            plot_sky_targets(outdir, spallfile, splog, nobs=True)

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
        cc = False
        if custom is not None:
            cc= True
        elif field is not None:
            field = field_to_string(field)
        if field is not None and mjd is not None:
            if not cc:
                field = field_to_string(field)
            mjd = str(mjd)
            specfull_dir =  field_spec_dir(indir, run2d,field, mjd, epoch=epoch,
                                           custom = cc, custom_name=custom)
            spallfile     = ptt.join(specfull_dir, 'spAll-'+field+'-'+mjd+'.fits.gz')
            spalllitefile = ptt.join(specfull_dir, 'spAll-lite-'+field+'-'+mjd+'.fits.gz')
            splinefile    = ptt.join(specfull_dir, 'spAllLine-'+field+'-'+mjd+'.fits.gz')
            spAlldatfile  = ptt.join(specfull_dir, 'spAll-'+field+'-'+mjd+'.dat.gz')
        else:
            fflags = []
            dflags = [indir]
            
            if run2d is not None:
                fflags.append(run2d)
                dflags.extend([run2d, 'summary'])
            else:
                dflags.append('summary')
            if custom is not None:
                dflags.append(custom)
                #spall_dir = ptt.join(indir, run2d, 'summary', fieldgroup(custom, custom=True))
                fflags.append(custom)
            elif epoch is True:
                dflags.append('epoch')
                fflags.append('epoch')
                #spAll_dir  = ptt.join(indir, run2d, 'summary','epoch')
            else:
                #spAll_dir  = ptt.join(indir, run2d)
                dflags.append('daily')
            
            spall_dir = ptt.join(*dflags)
            fflags = '-'.join(fflags)
            if len(fflags) > 0:
                fflags = '-'+fflags
            spallfile     = ptt.join(spall_dir, 'spAll'+fflags+'.fits.gz')
            spalllitefile = ptt.join(spall_dir, 'spAll-lite'+fflags+'.fits.gz')
            splinefile    = ptt.join(spall_dir, 'spAllLine'+fflags+'.fits.gz')
            spAlldatfile  = ptt.join(spall_dir, 'spAll'+fflags+'.dat.gz')

    if dev:
        spallfile = spallfile.replace('spAll','spAll_dev')
        spalllitefile = spalllitefile.replace('spAll','spAll_dev')
        splinefile = splinefile.replace('spAllLine','spAllLine_dev')
        spAlldatfile = spAlldatfile.replace('spAll','spAll_dev')
    
    return(spallfile, spalllitefile, splinefile, spAlldatfile)

def write_spAll(spAll, spline, spAll_lite, indir, run2d, datamodel, line_datamodel,
                epoch=False, dev=False, outroot=None, field=None, mjd = None,
                verbose=False, clobber=False, silent=False, custom = None,
                allsky = False, SDSSC2BV = '', tmpext = '.tmp', bkup = False):
 
    spallfile, spalllitefile, splinefile, spAlldatfile = build_fname(indir, run2d, outroot=outroot,
                                                                     field=field, mjd=mjd, dev=dev,
                                                                     epoch=epoch, custom=custom, allsky=allsky)

    drop_cols = None
    date = time.ctime()
    exists = ptt.exists(spallfile) if not clobber else False
    if spAll is None: return(spallfile)
    if bkup:
        if ptt.exists(spallfile):
            bkup_str = datetime.strptime(fits.getheader(spallfile,0)['DATE'],'%c').isoformat()
            
            if not ptt.exists(f"{spallfile}.bkup-{bkup_str}"):
                shutil.copy2(spallfile,f"{spallfile}.bkup-{bkup_str}")
                shutil.copy2(spalllitefile,f"{spalllitefile}.bkup-{bkup_str}")
                shutil.copy2(splinefile,f"{splinefile}.bkup-{bkup_str}")

    
    if not exists:
        hdul = merge_dm(ext = 'Primary', hdr = {'RUN2D':run2d,
                                                'Date':time.ctime(),
                                                'SDSSC2BV': SDSSC2BV},
                        dm = datamodel, splog=splog, verbose=verbose)
        splog.info('Formatting spAll table')
        spAll = merge_dm(table=spAll, ext = 'SPALL', name = 'SPALL', dm = datamodel,
                         splog=splog, drop_cols=drop_cols, verbose=verbose)
        makedirs(ptt.dirname(spallfile), exist_ok=True)
        splog.log('Writing '+spallfile)
        fits.HDUList([hdul,spAll]).writeto(spallfile.replace('.gz',tmpext+'.gz'),
                                           overwrite=True, checksum=True)
        hdul = None
        if spAll_lite is not None:
            spAll = Table(spAll.data)
            spAll.remove_columns(spAll_lite.colnames)#['ASSIGNED','ON_TARGET','VALID','DECOLLIDED',
                                 #         'MOON_DIST','MOON_PHASE','CARTON_TO_TARGET_PK'])
            spAll = hstack([spAll, spAll_lite], join_type = 'exact')
        else:
            spAll = None
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
            splog.info('Formatting spAll-lite table')
            spAll = merge_dm(table=spAll, ext = 'SPALL_lite', name = 'SPALL', dm = datamodel,
                                  splog=splog, drop_cols=drop_cols, verbose=verbose)
            splog.log('Writing '+spalllitefile)
            fits.HDUList([hdul_lite, spAll]).writeto(spalllitefile.replace('.gz',tmpext+'.gz'),
                                                          overwrite=True, checksum=True)
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
            splog.info('Formatting spAllLine table')
            spline = merge_dm(table=spline, ext = 'spZline', name = 'SPLINE', dm = line_datamodel,
                              splog=splog, drop_cols=drop_cols, verbose=verbose)
            splog.log('Writing '+splinefile)
            fits.HDUList([hdul_line,spline]).writeto(splinefile.replace('.gz',tmpext+'.gz'),
                                                     overwrite=True, checksum=True)
            del hdul_line, spline
            gc.collect()
        elif not silent:
            splog.log('Skipping '+splinefile+' (exists)')

    try:
        rename(spallfile.replace('.gz',tmpext+'.gz'),spallfile)
    except:
        pass
    try:
        rename(spalllitefile.replace('.gz',tmpext+'.gz'),spalllitefile)
    except:
        pass
    try:
        rename(splinefile.replace('.gz',tmpext+'.gz'),splinefile)
    except:
        pass

    #spAll.write(spAlldatfile, format='ascii.fixed_width_two_line')
    
    return(spallfile)
    
