#!/usr/bin/env python3
from boss_drp import idlspec2d_dir, MOUNTAIN, database_profile
from boss_drp.field.generations import generations
from boss_drp.utils.splog import splog, splog_name
from boss_drp.field import field_to_string
from boss_drp.utils import match as wwhere
from boss_drp.utils import (merge_dm, load_env, HiddenPrints)
from boss_drp.prep.readfibermaps.db_tools import *
try:
    from boss_drp.prep.GetconfSummary import find_confSummary, find_plPlugMapM
except:
    pass

from astropy.io import fits
from astropy.table import Table, vstack, join, Column, MaskedColumn, unique
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
import astropy.units as u
from glob import glob
import os.path as ptt
import os
import time
import numpy as np
import warnings
from time import sleep
from pydl.pydlutils.yanny import read_table_yanny, yanny
from pydl.pydlutils import sdss
from pydl import uniq

if not MOUNTAIN:
    try:
        from dustmaps.bayestar import BayestarQuery
        from dustmaps.sfd import SFDQuery
        from dustmaps.edenhofer2023 import Edenhofer2023Query as E3D2023Query
        from boss_drp.prep.simple_dust_2023 import simple_dust_2023
    except:
        splog.info('ERROR: dustmaps is not installed')
    
    try:
        from sdssdb.peewee.sdss5db.targetdb import database
        import sdssdb
        splog.add_external_handlers(sdssdb.log.name)
        test = database.set_profile(database_profile)

        if not test:
            splog.info('WARNING: No SDSSDB access - Defaulting to no_db')
            no_db_poss = True
        else:
            SDSSDBVersion = os.getenv('SDSSDB_VER',None)
            if SDSSDBVersion is None:
                SDSSDBVersion=sdssdb.__version__
    except:
        splog.info('WARNING: No SDSSDB access - Defaulting to no_db')
        no_db_poss = True
    else:
        no_db_poss = False
    from sdss_semaphore.targeting import TargetingFlags
    try:
        from  sdss_semaphore.targeting import logger as sem_log
        splog.add_external_handlers(sem_log.name)
    except:
        pass
else:
    no_db_poss=False
        
pratio = np.asarray([2.085, 2.085, 2.116, 2.134, 2.135])
chunkdata = None

try:
    bitmask = sdss.set_maskbits(maskbits_file=ptt.join(os.getenv('IDLUTILS_DIR'),"data/sdss/sdssMaskbits.par")) # Always call set_maskbits to set the latest version
except:
    pass


def readfibermaps(spplan2d=None, topdir=None, clobber=False, SOS=False, no_db=False, fast=False,
                  datamodel = None, SOS_opts=None, release = 'sdsswork', remote = False,
                  logger=None, V_TARG='*'):
    args = {'spplan2d':spplan2d, 'topdir':topdir, 'clobber':clobber,
          'SOS':SOS, 'no_db':no_db, 'fast':fast, 'datamodel': datamodel,
          'SOS_opts':SOS_opts, 'release': release, 'remote': remote,
          'logger':logger, 'V_TARG':V_TARG}
    if no_db_poss:
        no_db = True
    no_remote = not remote
    legacy=plates=fps=False
    if SOS is False:
        if datamodel is None:
            datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spfibermap_dm.par')
        if not ptt.exists(spplan2d):
            splog.info('spplan2d file '+spplan2d+' does not exists')
            exit()

        if topdir is None:
            topdir = ptt.dirname(spplan2d)
        plan = read_table_yanny(spplan2d,'SPEXP')
        plan.convert_bytestring_to_unicode()
        obs = plan.meta['OBS']
        run2d = plan.meta['RUN2D']
        mjd = int(plan.meta['MJD'])
        field = plan.meta['fieldname']
        
        if int(field) < generations.get('legacy', field=True)[1]:
            legacy = True
        elif int(field) < generations.get('legacy',field=True)[1]:
            plates = True
        else:
            fps = True
        spFibermap = 'spfibermap-'+field_to_string(field)+'-'+str(mjd)+'.fits'

    else:
        if not ptt.exists(topdir): os.makedirs(topdir, exist_ok=True)
        if datamodel is None:
            datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spfibermap_SOS_dm.par')

        fibermap_files = [SOS_opts['confSummary']]

        try:
            fibermap = read_table_yanny(fibermap_files[0],'FIBERMAP')
        except:
            splog.info(f"KeyError: 'No table named FIBERMAP in {fibermap_files[0]}!'")
            if logger is  None:
                splog.close()
            return
        
        
        run2d = os.getenv('RUN2D')
        obs = fibermap.meta['observatory']
        field = fibermap.meta['field_id']
        if int(field) == -999:
            field = 0
        confid = fibermap.meta['configuration_id']
        fps = True
        if SOS_opts['mjd'] is None:
            mjd = fibermap.meta['MJD']
        else:
            mjd = SOS_opts['mjd']
        spplan2d = 'SOS'
        ccd = SOS_opts['ccd']
        
        plan = Table(names = ('confid', 'fieldid', 'mjd', 'mapname', 'exptime'), dtype=(object, object, int, object, object ))

        spFibermap = 'spfibermap-'+field_to_string(field)+'-'+str(mjd)+'-'+ccd+'.fits'
        if ptt.exists(spFibermap):
            meta = fits.getdata(spFibermap,1)
            for row in meta:
                plan.add_row((row['CONFIGURATION_ID'],    row['MJD'], int(row['MJD']), row['CONFIGURATION_ID'], None))
            meta = None
        plan.add_row(        (confid,    field, int(mjd), confid, None))
            
        if topdir is None:
            topdir = ptt.join('/data', 'boss', 'sos', str(mjd))
        

        
    if not SOS:
        splog.open(logfile = ptt.join(topdir, spFibermap.replace('.fits','.log')), append= (not clobber))
        splog.info('Log file '+ptt.join(topdir, spFibermap.replace('.fits','.log'))+' opened '+ time.ctime())
    else:
        if logger is not None:
            if SOS_opts['log']:
                splog.add_file(ptt.join(SOS_opts['log_dir'], spFibermap.replace('.fits','.log')))
        else:
            if SOS_opts['log']:
                splog.open(logfile = ptt.join(SOS_opts['log_dir'], spFibermap.replace('.fits','.log')), append= (not clobber))

        splog.info('readfibermaps started at '+ time.ctime())
        try:
            with fits.open(ptt.join(topdir, spFibermap), mode='update', checksum=True, output_verify="fix") as hdul:
                hdul.verify('fix')
                hdul.flush()
            clobber = False
        except:
            clobber = True
            
    if not clobber:
        if (ptt.exists(ptt.join(topdir, spFibermap))):
            try:
                summary = Table(fits.getdata(ptt.join(topdir, spFibermap),'SUMMARY'))
                for col in summary.colnames:
                    summary[col] = summary[col].astype(object)
            except:
                splog.info('Failure reading SUMMARY from '+ptt.join(topdir, spFibermap))
                summary = merge_dm(ext = 'Summary', name = 'Summary', hdr=None, table = None, dm = datamodel)
                summary = Table(summary.data)
        else:
            summary = merge_dm(ext = 'Summary', name = 'Summary', hdr=None, table = None, dm = datamodel)
            summary = Table(summary.data)
    else:
        summary = merge_dm(ext = 'Summary', name = 'Summary', hdr=None, table = None, dm = datamodel)
        summary = Table(summary.data)
    

    plans = []
    hdul = None
    for i, row in enumerate(plan):
        if row['mapname'] in plans:
            continue
        plans.append(row['mapname'])
        
        if SOS:
            fibermap_file = fibermap_files[0]
        else:
            if fps:
                fibermap_file = find_confSummary(row['mapname'], obs=obs, release=release)
                if fibermap_file is None:
                    continue
            else:
                fibermap_file = find_plPlugMapM(str(mjd), row['fieldid'], row['mapname'], release=release)
                if fibermap_file is None:
                    continue

     
        if clobber is False:
            #display(summary)
            if ptt.basename(fibermap_file) in summary['EXTNAME'].data:#.decode():
                splog.info(ptt.basename(fibermap_file)+' is previously run')
                continue
            else:
                splog.info('Running: '+ (fibermap_file))
        else:
            splog.info('Running: '+ (fibermap_file))

        
        fibermap, hdr = buildfibermap(fibermap_file, run2d, obs, field, mjd, exptime = row['exptime'],
                                      fast=fast, fps=fps, plates=plates, legacy=legacy, SOS=SOS,
                                      no_db=no_db, indir = ptt.dirname(fibermap_file), release=release,
                                      no_remote=no_remote, V_TARG=V_TARG)
        
        if hdul is None:
            if (clobber) or (not ptt.exists(ptt.join(topdir, spFibermap))):
                hdu = merge_dm(ext = 'Primary', hdr = {'FIELD':field,'MJD':mjd, 'OBS':obs,'SPPLAN2D':ptt.basename(spplan2d),
                                                        'FPS':fps, 'Plate':plates,'Legacy':legacy }, dm = datamodel)
                
                hdul = fits.HDUList([hdu])
            else:
                try:
                    hdul = []
                    with fits.open(ptt.join(topdir, spFibermap)) as hdul_in:
                        #hdul=hdul_in.copy()
                        for hdu in hdul_in:
                            hdul.append(hdu.copy())
                    hdul = fits.HDUList(hdul)
                except:
                    os.remove(ptt.join(topdir, spFibermap))
                    args['clobber'] = True
                    splog.info(f'Failure opening {ptt.join(topdir, spFibermap)}... Clobbering and rerunning')
                    readfibermaps(**args)

        new_row = {}
        for key in hdr.keys():
            skip=False
            if key.upper() in ['MINSTDINBLOCKBOSS_SHARED', 'PLATEDESIGNVERSION', 'NINPUTS', 'PRIORITY', 'EVILSCAN']:
                continue
            for lcol in ['GUIDER_COEFF', 'PLATEINPUT']:
                if lcol in key.upper(): 
                    skip =True 
                    continue
            if skip is False:
                new_row[key.upper()] = [str(hdr[key])]
        new_row['EXTNAME'] = [ptt.basename(fibermap_file)]
        new_row['PLUGDIR'] = [ptt.dirname(fibermap_file)]
        sav_sum = summary.copy() if i !=0 else None
            
        sav_sum = merge_dm(ext = 'Summary', name = 'Summary', hdr=None, table = Table(new_row), dm = datamodel, old_tab = sav_sum)
        summary = Table(sav_sum.data)

        #(Table(sav_sum.data))
        if i == 0:
            exts= []
            for hdu in hdul[1:]:
                exts.append(hdu.header['EXTNAME'])
            if 'SUMMARY' in exts:
                hdul['SUMMARY'] = sav_sum
            else:
                hdul.append(sav_sum)
        else:
            hdul['SUMMARY'] = sav_sum
        
        hdul.append(merge_dm(ext = 'X', name = ptt.basename(fibermap_file), hdr=None, table = fibermap, dm = datamodel))
        
        splog.info('------------------------------------------------------------')
    if hdul is not None:
        while True:
            try:
                hdul.writeto(ptt.join(topdir, spFibermap),overwrite=True, checksum=True)
                break
            except:
                sleep(1)
                
            
    splog.info('Successful completion of readfibermaps at '+ time.ctime())
                
    splog.info('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    if not SOS:
        splog.close()
    else:
        if logger is not None:
            splog.close_file()
        else:
            splog.close()

 
 
def buildfibermap(fibermap_file, run2d, obs, field, mjd, exptime=None, indir=None,
                  fps=False, plates=False, legacy=False, SOS=False, fast=False,
                  no_db=False, release='sdsswork', no_remote=False, V_TARG='*'):
    
    if no_db is True:
        splog.info('Reading '+fibermap_file+' without DB access')

    if plates or legacy:
        fibermap = read_table_yanny(fibermap_file, 'PLUGMAPOBJ')
        fibermap.convert_bytestring_to_unicode()

        hdr = fibermap.meta
        
        for col in fibermap.colnames:
            if not ((fibermap[col].dtype  in [float, int, np.dtype('int16'), np.dtype('float32'), bool])):
                fibermap[col] = fibermap[col].astype(object)
                
        fibermap = readPlateplugMap(fibermap_file, fibermap, mjd, SOS=SOS, fast=fast,
                     exptime=exptime, plates=plates, legacy=legacy, no_db=no_db, indir=indir,
                     release=release, no_remote=no_remote, V_TARG=V_TARG)
        fibermap['fiber_offset'] = 0

    elif fps:
        fibermap = read_table_yanny(fibermap_file, 'FIBERMAP')
        fibermap.convert_bytestring_to_unicode()

        hdr = fibermap.meta
        for col in fibermap.colnames:
            dcol = fibermap[col]
            if not ((fibermap[col].dtype  in [float, int, np.dtype('int16'), np.dtype('float32'), bool])):
                fibermap[col] = fibermap[col].astype(object)
            elif (fibermap[col].dtype in [float, np.dtype('float32')]):
                dcol[dcol.data == -999] = np.nan
        fibermap = calcWokOffset(fibermap, fibermap_file)
        fibermap = readFPSconfSummary(fibermap, mjd, sos=SOS, no_db = no_db, fast=fast,V_TARG=V_TARG,
                                      release=release, no_remote=no_remote)
        hdr = fibermap.meta
        fibermap = flag_offset_fibers(fibermap)
        fibermap['SCI_EXPTIME'] = np.nan
        hdr['CARTRIDGEID'] = 'FPS-S' if hdr['observatory'] == 'LCO' else 'FPS-N'
    if 'TOO' not in fibermap.colnames and 'too' not in fibermap.colnames:
        fibermap['TOO'] = 0
    fibermap.meta = {}
    
    return(fibermap, hdr)
    


def flag_offset_fibers(fibermap):
    splog.info('Flagging offset Fibers')
    if not bool(int(fibermap.meta['is_dithered'])):
        #offsets = np.zeros(len(fibermap), dtype=bool)
        #fibermap.add_column(offsets,name='fiber_offset')
        foff = fibermap['fiber_offset']

        indx = np.where(np.logical_or((fibermap['delta_ra'].data != 0), (fibermap['delta_dec'].data != 0)))[0]
        foff[indx] = 1
    return(fibermap)


def mags2Flux(fibermap, correction):
    splog.info('Converting Magnitudes to Fluxes')
    flag = np.ones(len(fibermap), dtype=bool)
   
    if 'CatDB_mag' in fibermap.colnames:
        mags = fibermap['CatDB_mag'].data
        optical_prov = fibermap['optical_prov'].data.astype(str)
    else:
        mags = fibermap['mag'].data
        optical_prov = np.full(len(fibermap), 'fiber2mag')
    
    calibflux = fibermap['calibflux']
    cflux =  fibermap['calibflux'].data
    cfluxivar =  fibermap['calibflux_ivar'].data
    
    factor = np.exp(-correction/2.5 * np.log(10)) #AB correction factor for flux
    pratio = [2.085, 2.085, 2.116, 2.134, 2.135]# ratio of fiber2flux
    splog.info('PSF/fiber flux ratios = '+str(pratio))

    psf = np.isin(optical_prov, [x for x in list(set(optical_prov)) if 'psf' in x.lower()])

    for ifilt in range(5):
        ibad = np.where((cflux[:,ifilt] == 0) & (mags[:,ifilt] > 0) &
                        (mags[:,ifilt] < 50) & (flag == 1))[0]
        if len(ibad)>0:
            splog.info('Using plug-map fluxes for '+str(len(ibad))+' values in filter '+str(ifilt))

        ibad = np.where((cflux[:,ifilt] == 0) & (mags[:,ifilt] > 0) &
                        (mags[:,ifilt] < 50) & (flag == 1) & (psf == False))[0]
        if len(ibad)>0:
            cflux[ibad,ifilt] = np.power(10.,((22.5 - mags[ibad, ifilt]) / 2.5))*pratio[ifilt]
            cfluxivar[ibad,ifilt] = 0
            
        ibad = np.where((cflux[:,ifilt] == 0) & (mags[:,ifilt] > 0) &
                        (mags[:,ifilt] < 50) & (flag == 1) & (psf))[0]

        if len(ibad)>0:
            cflux[ibad,ifilt] = np.power(10.,((22.5 - mags[ibad, ifilt]) / 2.5))
            cfluxivar[ibad,ifilt] = 0
            
            
        #------------
        # Apply AB corrections to the CALIBFLUX values (but not to MAG)

    fibermap['calibflux']      = cflux * factor
    fibermap['calibflux_ivar'] = cfluxivar / np.power(factor,2)

    return(fibermap)


def psf2Fiber_mag(fibermap, plates=False, legacy=False):
    splog.info('Calculating PSFmags and fiber2mags from CatalogDB magnitudes')
    if plates or legacy:
        fibermap.add_columns([['fiber2mag'],[np.full(5,np.nan)]], names = ['optical_prov', 'fiber2mag'])

    else:
        fibermap.add_columns([[np.full(5,np.nan)],[np.full(5,np.nan)],[np.full(5,np.nan)]], names = ['CatDB_mag', 'fiber2mag','PSFmag'])

    if plates:
        magcol = fibermap['mag']
        magd = magcol.data
        for col in range(magd.shape[1]):
            mf = magd[:,col]
            mf[np.where(mf == -10)[0]] = np.nan
            mf[np.where((mf <= -9.999990) & (mf > -10.00001))[0]] = np.nan
            mf[np.where((mf <= -900.))[0]] = np.nan
            mf[np.where((mf <= 0.))[0]] = np.nan
            magd[:,col] = mf
        macol = magd

    fibermap['CatDB_mag'] = fibermap['mag'].data.copy()
    fibermap['fiber2mag'] = fibermap['mag'].data.copy()
    fibermap['PSFmag']    = fibermap['mag'].data.copy()
 
    optical_prov = fibermap['optical_prov'].data.astype(str)
    psf_optical_prov = [x for x in list(set(optical_prov)) if 'psf' in x.lower()]
    fiber2_optical_prov = [x for x in list(set(optical_prov)) if 'fiber2mag' in x.lower()]

    magcol = fibermap['mag']
    PSFmag = fibermap['PSFmag']
    fiber2mag = fibermap['fiber2mag']

    imatch = np.where(np.isin(optical_prov, psf_optical_prov))[0]
    if len(imatch) > 0:
        mag = fibermap['mag'].data.copy()
        mag = mag+2.5*np.log10(pratio)
        magcol[imatch]  = mag[imatch]
        fiber2mag[imatch] = mag[imatch]

    imatch = np.where(np.isin(optical_prov, fiber2_optical_prov))[0]
    if len(imatch) > 0:
        mag = fibermap['mag'].data.copy()
        mag = mag-2.5*np.log10(pratio)
        PSFmag[imatch] = mag[imatch]

    imatch = np.where(~np.isin(optical_prov, psf_optical_prov+fiber2_optical_prov))
    if len(imatch) > 0:
        fiber2mag[imatch] = np.full(5,np.nan)
        PSFmag[imatch]    = np.full(5,np.nan)
        
    imatch = np.where(np.isin(optical_prov, [x for x in list(set(optical_prov)) if 'other' in x.lower()]))[0]
    if len(imatch) > 0:
        fiber2mag[imatch] = np.full(5,np.nan)
        PSFmag[imatch] = np.full(5,np.nan)

    imatch = np.where(optical_prov == '')[0]
    if len(imatch) > 0:
        fiber2mag[imatch] = np.full(5,np.nan)
        PSFmag[imatch] = np.full(5,np.nan)
        magcol[imatch]    = np.full(5,np.nan)

    return(fibermap)


def get_survey(fibermap, plates= False, legacy = False):
    splog.info('Getting Survey Names')
    test = 'FIRSTCARTON' if ((plates is True) or (legacy is True)) else 'program'
    
    if 'survey' not in fibermap.colnames:
        fibermap.add_column(Column('', dtype=object, name='survey'))

    test_col = fibermap[test].astype(str).data
    surv = fibermap['survey']
    for survey in ['MWM', 'BHM', 'SKY', 'OPS', 'OPEN_FIBER', 'COMMISSIONING']:
        imatch = np.where(np.isin(test_col, [x for x in list(set(test_col)) if survey.upper() in x.upper()]))[0]
        if survey == 'SKY':
            survey = 'OPS'
        if len(imatch) > 0:
            surv[imatch] = survey
    survey = fibermap['survey'].data.astype(str)
    fibermap.remove_column('survey')
    fibermap.add_column(survey, name='survey')
    return(fibermap)


def catid_from_ra_dec(ra_deg, dec_deg):
    ra_deg = np.asarray(ra_deg, dtype=float)
    dec_deg = np.asarray(dec_deg, dtype=float)

    # RA: degrees -> hours
    ra_hours = ra_deg / 15.0
    rah = np.floor(ra_hours).astype(int)
    ram = np.floor((ra_hours - rah) * 60.0).astype(int)
    ras = ((ra_hours - rah) * 60.0 - ram) * 60.0

    # Dec: keep sign separate
    dec_sign = np.where(dec_deg < 0, "-", "+")
    dec_abs = np.abs(dec_deg)
    ded = np.floor(dec_abs).astype(int)
    dem = np.floor((dec_abs - ded) * 60.0).astype(int)
    des = ((dec_abs - ded) * 60.0 - dem) * 60.0

    return np.array([
        f"u{h:02d}{m:02d}{s:04.1f}{sg}{d:02d}{mm:02d}{ss:04.1f}"
        for h, m, s, sg, d, mm, ss in zip(rah, ram, ras, dec_sign, ded, dem, des)
    ])


def NoCatid(fibermap, plates = False, legacy = False):
    splog.info('Building pseudo IDs for fibers without catalogIDs')
    if (not plates) and (not legacy):
        ra  = fibermap['racat'].data
        dec = fibermap['deccat'].data
    else:
        ra  = fibermap['ra'].data
        dec = fibermap['dec'].data
    
    badidx = np.where(np.isnan(ra) | (np.isnan(dec)))[0]
    if len(badidx) > 0:
        ra[badidx]  = fibermap[badidx]['ra'].data
        dec[badidx] = fibermap[badidx]['dec'].data
        badidx = np.where(np.isnan(ra) | (np.isnan(dec)))[0]
        if len(badidx) > 0:
            ra[badidx]  = 0.0
            dec[badidx] = 0.0

    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    dummy_catid = catid_from_ra_dec(ra, dec)

    for c in badidx:
        pairidx = np.where(fibermap['holeId'] == fibermap[c]['holeId'])[0]
        if fibermap[c]['fiberType'] == 'APOGEE':
            pairidx1 = np.where(fibermap[pairidx]['fiberType'] == 'BOSS')[0]
            pairidx = pairidx[pairidx1[0]]
        elif  fibermap[c]['fiberType'] == 'BOSS':
            pairidx1 = np.where(fibermap[pairidx]['fiberType'] == 'APOGEE')[0]
            pairidx = pairidx[pairidx1[0]]
        elif  fibermap[c]['fiberType'] == 'METROLOGY':
            pairidx1 = np.where(fibermap[pairidx]['fiberType'] == 'BOSS')[0]
            pairidx = pairidx[pairidx1[0]]
        if pairidx in badidx:
            dummy_catid[c] = 'u'+str(fibermap[c]['fiberId'])
        else:
            if (not plates) and (not legacy):
                sc = SkyCoord(fibermap[pairidx]['racat'], fibermap[pairidx]['deccat'], frame='icrs', unit='deg')
            else:
                sc = SkyCoord(fibermap[pairidx]['ra'], fibermap[pairidx]['dec'], frame='icrs', unit='deg')
            sc = sc.to_string('hmsdms', sep='',precision=1)
            dummy_catid[c] = 'u'+sc.replace(' ','')
        
    if plates or legacy:
        iunassigned = np.where(fibermap['icatalogid'] == 0)[0]
        if len(iunassigned) > 0:
            fibermap_boss=fibermap[iunassigned]
            boss_dummpycat=dummy_catid[iunassigned]
            fibermap_boss['catalogid'] = boss_dummpycat
            fibermap_boss['icatalogid'] = -999
            fibermap[iunassigned]=fibermap_boss
    else:  
        All_apogee_fibers = np.where(fibermap['fiberType'] == 'APOGEE')[0]
        iunassigned_boss = np.where((fibermap['fiberType'] == 'BOSS') & (fibermap['assigned'] == 0))[0]
        if len(iunassigned_boss) > 0:
            fibermap_boss=fibermap[iunassigned_boss]
            boss_dummpycat=dummy_catid[iunassigned_boss]
            fibermap_boss['catalogid'] = boss_dummpycat
            fibermap_boss['icatalogid'] = -999
            fibermap[iunassigned_boss]=fibermap_boss

    return(fibermap)

def calcOffset(fibermap, obs_epoch):
    splog.info('Calculating Fiber Offsets in Sky Frame')
    try:
        catCord = SkyCoord(ra=np.ma.filled(fibermap['racat'])*u.deg,
                           dec=np.ma.filled(fibermap['bb'].data)*u.deg,
                           pm_ra_cosdec=np.ma.filled(fibermap['pmra'].data)*u.mas/u.yr,
                           pm_dec = np.ma.filled(fibermap['pmdec'].data)*u.mas/u.yr,
                           obstime=Time(np.nan_to_num(fibermap['coord_epoch'].data, nan=2016.0).astype(float), format='jyear'))
        obsCord = SkyCoord(ra=np.ma.filled(fibermap['racat'])*u.deg,
                           dec=np.ma.filled(fibermap['bb'].data)*u.deg,
                           obstime=Time(float(obs_epoch), format='jd'))
                       
        total_offset = catCord.separation(obsCord).arcsec
        fibermap.add_column(total_offset, name = 'MeasuredOffset')
    except:
        splog.info('Failure Calculating Fiber Offsets... Assuming fibers on Target')
    return(fibermap)
    
def calcWokOffset(fibermap, fibermap_file):
    splog.info('Calculating Fiber Offsets in Wok Frame')
    if 'confSummaryF' in fibermap_file:
        fibermap_pre = read_table_yanny(fibermap_file.replace('confSummaryF','confSummary'), 'FIBERMAP')
        fibermap_pre.meta = {}
        fibermap_pre.convert_bytestring_to_unicode()
        try:
            fibermap_pre = fibermap_pre[['positionerId', 'fiberType','xwok','ywok','zwok']]
        except:
            fibermap_pre['xwok'] = np.nan
            fibermap_pre['ywok'] = np.nan
            fibermap_pre['zwok'] = np.nan
            fibermap_pre = fibermap_pre[['positionerId', 'fiberType','xwok','ywok','zwok']]
            fibermap['xwok'] = np.nan
            fibermap['ywok'] = np.nan
            fibermap['zwok'] = np.nan
            
        for col in fibermap_pre.colnames:
            dcol = fibermap_pre[col]
            if not ((fibermap_pre[col].dtype  in [float, int, np.dtype('int16'), np.dtype('float32'), bool])):
                fibermap_pre[col] = fibermap_pre[col].astype(object)
            elif (fibermap_pre[col].dtype in [float, np.dtype('float32')]):
                dcol[dcol.data == -999] = np.nan

        
        fibermap = join(fibermap, fibermap_pre, keys=['positionerId', 'fiberType'], join_type='left', table_names=['', '_pre'])


        for col in ['xwok','ywok','zwok']:
            fibermap.rename_column(col+'_', col)
            fibermap.rename_column(col+'__pre', col+'_pre')

            try:
                if fibermap.masked:
                    if np.any(fibermap[col].mask):
                        mask = fibermap[col].mask
                        fibermap[col][mask] = fibermap[col][mask]
            except:
                pass

        xwok_pre = fibermap['xwok'].data
        ywok_pre = fibermap['ywok'].data
        fibermap['WokOffset'] = np.sqrt((fibermap['xwok'] - fibermap['xwok_pre'])**2 + (fibermap['ywok'] - fibermap['ywok_pre'])**2)

        fibermap_pre= None

    else:
        try:
            fibermap['xwok_pre']=fibermap['xwok']
            fibermap['ywok_pre']=fibermap['ywok']
            fibermap['zwok_pre']=fibermap['zwok']
            fibermap.add_column(0, name = 'WokOffset')
        except:
            fibermap['xwok'] = np.nan
            fibermap['ywok'] = np.nan
            fibermap['zwok'] = np.nan
            fibermap['xwok_pre']=fibermap['xwok']
            fibermap['ywok_pre']=fibermap['ywok']
            fibermap['zwok_pre']=fibermap['zwok']
            fibermap.add_column(np.nan, name = 'WokOffset')
    return(fibermap)

def readFPSconfSummary(fibermap, mjd, sos=False, no_db = False, fibermask = None, fast=False,
                       release='sdsswork', no_remote=False, V_TARG='*'):
    
    # The correction vector is here --- adjust this as necessary.
    # These are the same numbers as in SDSSFLUX2AB in the photoop product.
    
    correction = np.asarray([-0.042, 0.036, 0.015, 0.013, -0.002])
    nfiber = len(fibermap)
   
    fibermap.rename_column('catalogid', 'icatalogid')
    fibermap.add_column(Column(fibermap['icatalogid'].astype(object), name = 'catalogid'))

    fibermap = NoCatid(fibermap, plates = False, legacy = False)
    
    if fibermask is None:
        fibermask = np.zeros(nfiber, dtype=int)
    if len(fibermask) != nfiber:
        splog.info('Number of elements in FIBERMASK do not match NFIBER')
        exit()
    
    fibermap.add_columns([0,1,False,0,fibermap['category'].astype(object)], names = ['badstdmask', 'offsetid','fiber_offset','fibermask','objtype'])

    if ((bool(int(fibermap.meta['is_dithered']))) or (int(fibermap.meta['parent_configuration']) != -999)):
        iAssigned = np.where((fibermap['assigned'] == 1) & (fibermap['valid'] == 1))[0]
    else:
        iAssigned = np.where((fibermap['assigned'] == 1) & (fibermap['on_target'] == 1) & (fibermap['valid'] == 1))[0]


    fibermask = np.bitwise_or(fibermask.astype(int), sdss.sdss_flagval('SPPIXMASK', 'NOPLUG'))
    if len(iAssigned) > 0: 
        fibermask[iAssigned] = fibermask[iAssigned] - sdss.sdss_flagval('SPPIXMASK', 'NOPLUG')
    fibermap['fibermask']=fibermask

    if sos is True:
        fibermap = psf2Fiber_mag(fibermap, plates=False, legacy=False)
        hdr = fibermap.meta
        fibermap.add_column(Column(int(hdr['configuration_id']), name = 'configuration_id'))
        fibermap.add_column(Column(hdr['robostrategy_run'], name = 'targeting_vers', dtype=object))
        fibermap.add_column(Column(int(hdr['field_id']), name = 'fieldid'))
        fibermap.add_column(Column(int(hdr['MJD']), name = 'MJD'))
        fibermap.add_column(Column(float(hdr['raCen']), name = 'rafield'))
        fibermap.add_column(Column(float(hdr['decCen']), name = 'decfield'))
        fibermap.add_column(Column(float(0), name = 'redden_med'))
        fibermap.add_column(Column([np.zeros(3,dtype=float)], name = 'fibersn'))
        fibermap.add_column(Column([np.zeros(3,dtype=float)], name = 'synthmag'))
        fibermap.add_column(Column(float(0), name = 'hrmed'))
        fibermap.add_column(Column([np.zeros(5,dtype=float)], name = 'calibflux'))
        fibermap.add_column(Column([np.zeros(5,dtype=float)], name = 'calibflux_ivar'))
        
    mjd=int(fibermap.meta['MJD'])
    lco = True if fibermap.meta['observatory'].upper() == 'LCO' else False
    # Read calibObj or photoPlate photometry data
    if not sos:
        fieldid = fibermap.meta['field_id']
        ra_field = float(fibermap.meta['raCen'])
        dec_field = float(fibermap.meta['decCen'])
        obs_epoch = float(fibermap.meta['epoch'])
        
        fibermap=calibrobj(fibermap, fieldid, ra_field, dec_field, fps=True, 
                           lco=lco, design_id=fibermap.meta['design_id'], fast=fast,
                           RS_plan=fibermap.meta['robostrategy_run'], no_db=no_db,
                           release=release, no_remote=no_remote, V_TARG=V_TARG)
        fibermap = calcOffset(fibermap, obs_epoch)
    else:
        fibermap=mags2Flux(fibermap, correction)
    objtype = fibermap['objtype']
    
    program = fibermap['program'].data
    program = np.char.lower(program.astype(str))
    
    objtype[np.where(fibermap['objtype'].data == 'sky_boss')[0]] = 'SKY'
    objtype[np.where(fibermap['objtype'].data == 'ops_sky')[0]] = 'SKY'
    objtype[np.where(wwhere(program, '*ops_std*'))[0]] = 'SPECTROPHOTO_STD'


    fibermap = get_survey(fibermap)
    fibermap = fps_fibermapsort(fibermap)
    for col in ['carton', 'program_db', 'll', 'bb', 'rr', 'stdflag']:
        if col in fibermap.colnames:
            fibermap.remove_column(col)
    
    return(fibermap)
 
 
def calibrobj(fibermap, fieldid, rafield, decfield, design_id=None, 
              plates=False, legacy=False, fps=False, sos=False, lco=False, RS_plan=None, 
              no_db=False, mjd=None, indir=None, fast =False, release='sdsswork',
              no_remote=False, designmode=None,V_TARG='*'):

    # The correction vector is here --- adjust this as necessary.
    # These are the same numbers as in SDSSFLUX2AB in the photoop product.
    correction = np.asarray([-0.042, 0.036, 0.015, 0.013, -0.002])

    splog.info('Adding fields from calibObj file')
        
    hdr = fibermap.meta
    fibermap.add_column([np.full(4,np.nan)], name = 'WISE_MAG')
    fibermap.add_column([np.full(3,np.nan)], name = 'TWOMASS_MAG')
    fibermap.add_column([np.full(2,np.nan)], name = 'GUVCAT_MAG')
    if fps:
        fibermap.add_column(fibermap['firstcarton'].data, name = 'CARTONNAME')
    else:
        fibermap.add_column(fibermap['FIRSTCARTON'].data, name = 'CARTONNAME')
    fibermap= get_supplements(fibermap, designID=design_id, rs_plan = RS_plan,
                              fps= fps, fast=fast, release=release, designmode=designmode,
                              no_remote=no_remote, db = (not no_db), V_TARG=V_TARG)
    fibermap.add_column([np.zeros(5,dtype=float)], name = 'calibflux')
    fibermap.add_column([np.zeros(5,dtype=float)], name = 'calibflux_ivar')
    fibermap.add_column([np.zeros(5,dtype=int)], name = 'calib_status')
    #----------
    # Attempt to read the calibObj photometry data
    if legacy:
        tsobj = plug2tsobj(fieldid, plates=plates, legacy=legacy,  mjd=mjd, indir=indir)
        # Do not use the calibObj structure if more than 20% of the non-sky
        # objects do not have fluxes.
        if tsobj is not None:
            qexist = tsobj['psfflux'][2] != 0
            qsky = np.isin((fibermap['objtype'].data, [x for x in list(set(fibermap['objtype'].data)) if 'sky' in x.lower()]))
            splog.info('Matched '+str(int(sum((qsky == 0) & qexist)))+' of '+str(int(sum((qsky == 0))))+' non-SKY objects')
            if (sum((qsky == 0) & qexist) < 0.80*sum(qsky == 0)):
                splog, 'Discarding calibObj structure because < 80% matches'
                tsobj = None
        if tsobj is not None:
            # Propagate CALIB_STATUS information:
            if 'CALIB_STATUS' in tsobj.colnames:
                fibermap['calib_status'] = tsobj['calib_status']
            #Assume that all objects not called a 'GALAXY' are stellar objects
            qstar = not np.isin((fibermap['objtype'].data, [x for x in list(set(fibermap['objtype'].data)) if 'galaxy' in x.lower()]))
            istar = np.where(qstar & qexist)[0]
            nstar = len(istar)
            igal  = np.where((not qstar) & qexist)[0]
            if 'fiber2flux' in tsobj.colnames:
                pratio = [2.085, 2.085, 2.116, 2.134, 2.135]
            else:
                pratio = [1.343, 1.336, 1.354, 1.363, 1.367]
            if (nstar > 0):
                fibermap[istar]['calibflux'] = tsobj[istar]['psfflux']
                fibermap[istar]['calibflux_ivar'] = tsobj[istar]['psfflux_ivar']
            splog.info('PSF/fiber flux ratios = '+str(pratio))
            if (len(istar) > 0):
                    fibermap[igal]['calibflux']      = tsobj['fiberflux'][igal] * pratio
                    fibermap[igal]['calibflux_ivar'] = tsobj['fiberflux_ivar'][igal] / ((pratio)*(pratio))

            # Reject any fluxes based upon suspect PHOTO measurements, as indicated by the PHOTO flags.
            badbits2 = sdss.sdss_flagval('OBJECT2', 'SATUR_CENTER') | sdss.sdss_flagval('OBJECT2', 'INTERP_CENTER') | sdss.sdss_flagval('OBJECT2', 'PSF_FLUX_INTERP') 

            qgoodphot = (tsobj['flags2'] & badbits2) == 0
            fibermap['calibflux'] = fibermap['calibflux'] * qgoodphot
            fibermap['calibflux_ivar'] = fibermap['calibflux_ivar'] * qgoodphot
        else:
            splog.info('WARNING: No calibObj structure found for plate '+fieldid)

    if fps:
        fibermap = psf2Fiber_mag(fibermap,plates=False, legacy=False)
    fibermap = mags2Flux(fibermap, correction)
     
    return(fibermap)


def plug2tsobj(plateid, ra=None, dec=None, mjd=None, indir=None, dmin=2.0, 
               silent=False, plugmap = None):

    if not((type(plateid) == str) or (type(plateid) == int)):
        splog.info('PLATEID must be a scalar!')
        exit()
        
    if (mjd is not None):
        if not ((type(plateid) == str) or (type(mjd) == int)):
            splog.info('Number of elements in PLATEID and MJD must agree!')
            exit()

    platestr = field_to_string(plateid)

    if indir is None:
        indir = ptt.join(os.getenv('BOSS_SPECTRO_REDUX'),os.getenv('RUN2D'),platestr)
    
    if ra is not None:
        if len(ra) != len(dec):
            splog.info('Number of elements in RA and DEC must agree!')
            exit()

    #----------
    # First look for photoPosField files (if MJD set)
    filename = None
    if mjd is not None:
        mjdstr = str(int(mjd))
        shortname = 'photoPosPlate-'+platestr+'-'+mjdstr+'.fits'
        # Look in the output RUN2D directory first, then any subdirectories if not found
        filename = ptt.join(indir,shortname)
        qsorted = True
        if not ptt.exists(filename):
            filename = None

    #----------
    # Next look for calibPlateP file
    if filename is None:
        filename = 'calibPlateP-' + platestr + '.fits'
        filename = ptt.join(indir, filename)
        if ptt.exists(filename):
            qsorted = False
        else:
            filename = None

    #----------
    # Read the file
    if filename is None:
        splog.info( 'WARNING: photoPosPlate file not found for plate ' + platestr)
        return(None)
    splog, 'Reading object calibration file: ' + filename

    # Make certain that the file exists and is valid
    message = 0
    try: 
        tstemp = Table(fits.getdata(filename,1))
    except:
        splog.info('WARNING: calibObj file is empty: ' + filename)
        return(None)

    #----------
    # Sort the file if necessary

    if (qsorted): 
        return(tstemp)
    
    test_cols = []
    for column in tstemp.columns:
        test_cols.append([(column,tstemp[column].dtype)][0])
    tsobj = Table(data=np.zeros(len(ra), dtype=test_cols))

    #----------
    # Match objects by positions on the sky
    for iplug, r in enumerate(ra):
        # Assume that this object is non-existent if RA=DEC=0
        if ((ra[iplug] != 0 & dec[iplug] != 0) or not(np.isnan(ra[iplug]) or np.isnan(dec[iplug]))):
            c1 = SkyCoord(tstemp.matchra*u.deg, tstemp.matchdec*u.deg)
            c2 = SkyCoord(ra[iplug]*u.deg, dec[iplug]*u.deg)
            sep = c1.separation(c2).arcsec.value
            imin = np.argmin(sep)
            if min(sep) > dmin:
                if plugmap is not None:
                    if 'SKY' not in plugmap[iplug]['objtype']:
                         splog.info('Warning: Unmatched OBJTYPE='+ plugmap[iplug]['objtype'] + ' at RA='+ str(ra[iplug])+ ',DEC='+str(dec[iplug]))
        
                else: 
                    splog.info('Warning: Unmatched ' + ' at RA='+ str(ra[iplug])+ ',DEC='+str(dec[iplug]))
            else:
                tsobj[iplug] = tstemp[imin]
    return(tsobj)


def fps_fibermapsort(fibermap, add =True):
    splog.info('Sorting FPS Fibermap')
    if add:
        fibermap.add_column(fibermap['fiberId'].data, name = 'ConfFiberid')
    else:
        fibermap.add_column(fibermap['fiberId'].data, name = 'ConfFiberid_1')

    fibermap.sort(['fiberId'])
   
    fid = fibermap['fiberId']
    fibermap=vstack([fibermap[fibermap['fiberType'].astype(str) == 'BOSS'],
                     fibermap[fibermap['fiberType'].astype(str) == 'APOGEE'],
                     fibermap[fibermap['fiberType'].astype(str) == 'METROLOGY']])
    fid = fibermap['fiberId']
    fid[fibermap['fiberType'].data == 'BOSS'] = np.arange(1,len(fibermap[fibermap['fiberType'].data == 'BOSS'])+1)
    fibermap['fiberId'] = fid
    return(fibermap)


def plate_fibermapsort(fibermap, fibermask=None, plates = False):
    splog.info('Sorting Plate Fibermap')
    fibermap_raw = fibermap.copy()
    fibermap = fibermap[fibermap['holeType'] == 'OBJECT']
    
    if fibermask is None:
        fibermask = np.zeros(len(fibermap))
    elif len(fibermask) != len(fibermap):
        splog.info('Number of elements in FIBERMASK do not match NFIBER')
        exit()    
    
    badstdmask = np.zeros(len(fibermap), dtype=bool)
    fibermap.add_column(Column(badstdmask, name = 'badstdmask'))
    fibermap.add_column(Column(np.full(len(fibermap),2000.0), name='COORD_EPOCH'))
    if plates:
        programname = fibermap.meta['programname']
        fibermap.add_column(Column(fibermap['fiberId'].data, name = 'org_fiberid'))
        fibermap.add_column(Column(programname, name = 'program', dtype=object))
        stdmask = fibermap['badstdmask']
        if 'MWM' in programname.upper():
            for i, row in enumerate(fibermap):
                if row['objType'] == 'SPECTROPHOTO_STD':
                    if (row['mag'])[3] >= 18.0:
                        badstdmask[i] = 1
                        stdmask[i] = 1
        
        igood =      np.where((fibermap['fiberId'].data > 0) & (fibermap['spectrographId'].data == 1) & (badstdmask == 0))[0]
        iplugged =   np.where((fibermap['fiberId'].data > 0) & (fibermap['spectrographId'].data == 1))[0]
        igoodapoge = np.where((fibermap['fiberId'].data > 0) & (fibermap['spectrographId'].data == 2) & (badstdmask == 0))[0]
    else:
        igood =      np.where((fibermap['holeType'] == 'OBJECT') & (fibermap['fiberid'].data > 0))[0]
        iplugged =   igood
    
    if len(igood) == 0:
        splog.info('No fibers found in plugmap!')
        exit()

    # Set the appropriate fibermask bit if a fiber not found in plugmap file.
    # Do this by first setting all bits to 1, then unsetting the good ones.
    fibermask = np.bitwise_or(fibermask.astype(int), sdss.sdss_flagval('SPPIXMASK', 'NOPLUG'))
    #fibermaks = fibermask | sdss.sdss_flagval('SPPIXMASK', 'NOPLUG')
    fibermask[iplugged] = fibermask[iplugged] - sdss.sdss_flagval('SPPIXMASK', 'NOPLUG')
    fibermap.add_column(Column(fibermask, name ='fibermask'))
    fibermap.add_column(Column(False, name = 'fiber_offset'))
    fibermap.add_column(Column(0, name = 'MeasuredOffset'))
    fibermap.add_column(Column(0.0, name = 'delta_ra'))
    fibermap.add_column(Column(0.0, name = 'delta_dec'))
    
    

    # Fill in unplugged fibers with arbitrary entries, and assign
    # them a FIBERID.  After this, plugsort.fiberid should run from 1...nfiber
    sid = fibermap['spectrographId'].data
    fid = fibermap['fiberId'].data
    if plates:
        SP1  = fibermap[(sid == 1) & (fid > 0)]
        SP2  = fibermap[(sid == 2) & (fid > 0)]
        unassigned = fibermap[np.logical_or((fid <= 0), ((sid != 1) & (sid != 2)))]
    if plates:
        SP2['fiberId'] = SP2['fiberId'].data + 500
    
    missing = np.setdiff1d(np.arange(1,len(fibermap)+1),(vstack([SP1,SP2]))['fiberId'].data)        
    unassigned['fiberId'] = missing.astype(object)
    
    fibermap = vstack([SP1,SP2,unassigned])        
    fibermap.sort('fiberId')
    
    if plates:
        missing = np.setdiff1d(np.arange(1,500+1),SP1['fiberId'].data)
    splog.info('Number of missing fibers: '+str(len(missing)))
    splog.info('Number of Invalid Standards: '+str(sum(badstdmask)))

    return(fibermap)


def readPlateplugMap(plugfile, fibermap, mjd, SOS=False, V_TARG='*',
                     exptime=None, fibermask=None, plates=False,
                     legacy=False, no_db=False,indir=None, fast=False,
                     release='sdsswork', no_remote=False):

    # The correction vector is here --- adjust this as necessary.
    # These are the same numbers as in SDSSFLUX2AB in the photoop product.
    correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

    #----------
    # Trim to object fibers only, sort them, and trim to spectrographid

    fibermap = plate_fibermapsort(fibermap, fibermask=fibermask, plates=plates)
    #----------
    # Add the tags OFFSETID and SCI_EXPTIME for
    fibermap.add_columns([0,0.0], names = ['offsetid','SCI_EXPTIME'])

    offsetid    = fibermap['offsetid']
    SCI_EXPTIME = fibermap['SCI_EXPTIME']
    try:
        plugpoint = read_table_yanny(fibermap,'PLUGMAPPOINT')
        splog.info('Using OFFSETID and SCI_EXPTIME from PLUGMAPPOINT structure')
        for row in plugpoint:
            k = np.where((np.abs(fibermap['xFocal'].data - row['xFocal'].data) < 0.0001) &
                         (np.abs(fibermap['xFocal'].data - row['xFocal'].data) < 0.0001))[0]
            if len(k) > 0:
                offsetid[k[0]]    = row['offsetid']
                SCI_EXPTIME[k[0]] = row['SCI_EXPTIME']
    except:
        # Use default values
        fibermap['offsetid'] = 1


    fibermap['RACAT']  = fibermap['ra'].data
    fibermap['DECCAT'] = fibermap['dec'].data
    if (exptime is not None):
        iuniq = uniq(fibermap['offsetid'].data, np.sort(fibermap['offsetid'].data))
        exptot = np.sum(SCI_EXPTIME[iuniq].data)
        if (exptot > 0):
            splog.info('Rescaling SCI_EXPTIME values by '+str(exptime/exptot))
            SCI_EXPTIME = fibermap['SCI_EXPTIME'].data * exptime/exptot

    plateid = fibermap.meta['plateId']
    redden_med = np.asarray(fibermap.meta['reddeningMed'].split()).astype(float)
    if len(redden_med) != 5:
        splog.info('WARNING: Wrong number of elements for reddeningMed')
        redden_med = np.zeros(5,dtype=float)

    #----------
    # Append some information from the plateHoles file

    platelist_dir = os.getenv('PLATELIST_DIR')
    if platelist_dir is None:
        splog.info('ERROR: PLATELIST_DIR environmental variable must be set')
        exit()
    platefile = 'plateHoles-' + str(int(plateid)).zfill(6) + '.par'
    if platelist_dir is not None:
        thisfile = glob(ptt.join(platelist_dir, 'plates','*','*', platefile))
        if len(thisfile) > 0:
            plateholes = read_table_yanny(thisfile[0],'STRUCT1')
            plateholes.convert_bytestring_to_unicode()
            iobj = np.where(np.isin(plateholes['holetype'].data, [x for x in list(set(plateholes['holetype'].data)) if 'boss' in x.lower()]))[0]
            isort = np.full(len(fibermap), -1, dtype=int)
            for i, row in enumerate(fibermap):
                try:
                    isort[i] = np.where((plateholes['xfocal'] == row['xFocal']) & (plateholes['yfocal'] == row['yFocal']))[0]          
                except:
                    isort[i] = 0
            plateholes = plateholes[isort]            
            test_cols = []
            for column in plateholes.columns:
                test_cols.append([(column,plateholes[column].dtype)][0])
            blankhole = Table(data=np.zeros(1, dtype=test_cols))
            
            ibad = np.where(iobj == -1)
            for i in ibad:
                plateholes[i] = blankhole
            if plates:
                htags = ['SOURCETYPE','LAMBDA_EFF','ZOFFSET','BLUEFIBER', 
                         'BOSS_TARGET*','ANCILLARY_TARGET*', 'EBOSS_TARGET*', 
                         'CATALOGID','SDSSV_BOSS_TARGET*','FIRSTCARTON', 
                         'RUN','RERUN','CAMCOL','FIELD','ID', 'THING_ID_TARGETING']                
                
                fibermap.add_columns([0.0,0.0,0.0], names = ['Gaia_G_mag','BP_mag','RP_mag'])
                fibermap['Gaia_G_mag'] = plateholes['gaia_g']
                fibermap['BP_mag'] = plateholes['gaia_bp']
                fibermap['RP_mag'] = plateholes['gaia_rp']
            else:
                htags = ['SOURCETYPE','LAMBDA_EFF','ZOFFSET','BLUEFIBER', 
                         'BOSS_TARGET*','ANCILLARY_TARGET*', 'EBOSS_TARGET*', 
                         'RUN','RERUN','CAMCOL','FIELD','ID', 'THING_ID_TARGETING']
                
            for tag in list(htags):
                if '*' in tag:
                    htags.remove(tag)
                    n_htags = [x for x in list(set(plateholes.colnames)) if tag.replace('*','').lower() in x.lower()]
                    htags.extend(n_htags)
            htags = list(set(htags))
            for col in htags:
                fibermap.add_column(Column(plateholes[col.lower()], name = col))

            #- We never used washers < 175 microns
            ZOFFSET = fibermap['ZOFFSET']
            ii = np.where(fibermap['ZOFFSET'].data < 175)[0]
            ZOFFSET[ii] = 0

            #- Check opfiles/washers.par for overrides to ZOFFSET
            washers_file = ptt.join(idlspec2d_dir, 'opfiles','washers.par')
            washers = read_table_yanny(washers_file,'WASHERSTATUS')

            # extract fibermap file name to match header keyword NAME
            # plPlugMapM-5317-56000-01.par -> 5317-56000-01
            tmp = (ptt.basename(plugfile)).split('.')[0].split('-')
            plugname = '-'.join(tmp[1:])
            mjd = int(tmp[2])

            ii = np.where(washers['plugname'].data == plugname)[0]
            n = len(ii)
            if (n > 1):
                splog.info("ERROR: multiple washers.par entries for " + plugname)
                exit()
            if (n == 1):
                status = washers[ii[0]]['status']
                splog.info("INFO: washer ZOFFSET override "+ plugname+ " "+ status)
                if (status != 'Y'):
                    if (status == 'N'): 
                        ZOFFSET = 0.0
                    if (status == 'L'): 
                        ZOFFSET = (ZOFFSET.data == 0) * 300.0
                    if (status == 'T'): 
                        ZOFFSET = (ZOFFSET.data == 300) * 300.0
                    if (status == 'X'): 
                        splog.info("WARNING: We know that we don't know ZOFFSET washer status for "+ plugname)
                        splog.info("WARNING: setting washer ZOFFSET to default 0.0 for "+ plugname)
                        ZOFFSET = 0.0
            else:
                # No explicit override; check mjd before washers were available
                # don't print info for current plates to keep SoS quiet
                if (mjd < 56200): 
                    splog.info("INFO: no washers.par entry for "+ plugname)
                if (mjd < 55442):
                    splog.info("INFO: setting ZOFFSET=0 for MJD "+str(mjd)+ " < 55442")
                    ZOFFSET = 0.0
                # No more washers after plate 7184
                if (int(plateid) > 7184):
                    ZOFFSET = 0.0
                    
    if plates:
        #Check if plate is in list of plates to be corrected
        catfile_dir = ptt.join(idlspec2d_dir,'catfiles')
        if not ptt.exists(catfile_dir):
            splog.info('WARNING: No Catalogid Correction Files found in '+catfile_dir)
        else:
            catfile = glob(ptt.join(catfile_dir, 'Corrected_values_plate'+plateid.strip()+'_design*.fits'))
            if len(catfile) > 0:
                splog.info('Correcting Catalogid, Carton, and SDSS Magnitudes for '+plateid)
                catdata = Table(fits.getdata(catfile[0],1))
                catdata.convert_bytestring_to_unicode()
                
                fibermap.add_column(Column(0, name = 'SDSSV_APOGEE_TARGET0'))
                fibermap.add_column(Column(0, name = 'GRI_GAIA_TRANSFORM'))
                fibermap.add_column(Column(fibermap['CATALOGID'].data, name = 'org_catid'))
                
                
                catid = fibermap['CATALOGID']
                mag   = fibermap['mag']
                fcart = fibermap['FIRSTCARTON']
                at0   = fibermap['SDSSV_APOGEE_TARGET0']
                bt0   = fibermap['sdssv_boss_target0']
                gg0   = fibermap['GRI_GAIA_TRANSFORM']
                for j, row in enumerate(fibermap):
                    c1 = SkyCoord(catdata['RA'].data*u.deg, catdata['Dec'].data*u.deg)
                    c2 = SkyCoord(row['ra']*u.deg, row['dec']*u.deg)
                    sep = c1.separation(c2).arcsec
                    imin = np.argmin(sep)                    
                    
                    if sep[imin] > 2.0:
                        splog.info("WARNING: Catalogid Missmatch "+str(j+1)+" @RA DEC:"+str(row['ra'])+" "+str(row['dec']))
                        continue
                    else:
                        splog.info("Correcting Fiber "+str(j+1)+" Info @RA DEC:"+str(row['ra'])+" "+str(row['dec']))
                        catid[j] = catdata['Final_CatalogID'].data[imin]
                        mag[j,1] = catdata['gmag'].data[imin]
                        mag[j,2] = catdata['rmag'].data[imin]
                        mag[j,3] = catdata['imag'].data[imin]
                        mag[j,4] = catdata['zmag'].data[imin]
                        fcart[j] = catdata['First_Carton'].data[imin]
                        at0[j]   = catdata['APOGEE_Flag'].data[imin]
                        bt0[j]   = catdata['BOSS_Flag'].data[imin]
                        gg0[j]   = catdata['Transformation_Flag'].data[imin]
        fibermap.add_column(Column(-999, name='ASSIGNED'))
        fibermap.add_column(Column(-999, name='ON_TARGET'))
        fibermap.add_column(Column(-999, name='VALID'))
        fibermap.add_column(Column(-999, name='DECOLLIDED'))

    mag   = fibermap['mag']
    for i, row in enumerate(fibermap):
        indx = np.where(row['mag'].data == 0)[0]
        if len(indx) > 0:
            mag[np.where(row['mag'].data == 0)[0],i] = np.nan

    for col in  ['Gaia_G_mag','BP_mag','RP_mag']:
        mag   = fibermap[col]

        mag[np.where(fibermap[col].data == 0)[0]] = np.nan

    fibermap = psf2Fiber_mag(fibermap, plates=plates, legacy=legacy)
    
    fibermap.rename_column('CATALOGID', 'icatalogid')
    fibermap.add_column(fibermap['icatalogid'].astype(object), name = 'catalogid')

    fibermap = NoCatid(fibermap, plates = plates, legacy=legacy)

    #----------
    # Optionally add tags for SOS

    if SOS:
        fibermap.add_column(Column(fibermap.meta['cartridgeId'], name = 'cartid'))
        fibermap.add_column(Column(int(plateid), name = 'plateid'))
        fibermap.add_column(Column(int(fibermap.meta['tileId']), name = 'tileid'))
        fibermap.add_column(Column(float(fibermap.meta['raCen']), name = 'raplate'))
        fibermap.add_column(Column(float(fibermap.meta['decCen']), name = 'decplate'))
        fibermap.add_column(Column(float(redden_med), name = 'redden_med'))
        fibermap.add_column(Column(np.zeros(3, dtype=float), name = 'fibersn'))
        fibermap.add_column(Column(np.zeros(3, dtype=float), name = 'synthmag'))

    mag   = fibermap['mag']
    if plates:
        programname = fibermap.meta['programname']
        if 'eFEDS' in programname:
            psffibercor = [0.7978, 0.8138, 0.8230, 0.8235]
            spht = np.isin(fibermap['FIRSTCARTON'].data, [x for x in list(set(fibermap['FIRSTCARTON'].data)) if 'bhm_spiders_clusters-efeds' in x.lower()])
            for i, row in enumerate(fibermap):
                if 'bhm_spiders_clusters-efeds' in row['FIRSTCARTON']:
                    mag.data[i,1:] = mag.data[i,1:] - psffibercor
        fibermap.add_column(Column('Plates', name = 'cadence'))

        p2d = {'AQMES-Bonus':'dark_monit_plates',
               'AQMES-Medium':'dark_monit_plates',
               'AQMES-Wide':'dark_monit_plates',
               'eFEDS1':'dark_faint_plates',
               'eFEDS2':'dark_faint_plates',
               'eFEDS3':'dark_faint_plates',
               'MWM':'bright_time_plates',
               'MWM2':'bright_time_plates',
               'MWM2_sky':'bright_time_plates',
               'MWM3':'bright_time_plates',
               'MWM_30min':'bright_time_plates',
               'MWM_30min2':'bright_time_plates',
               'MWM_30min3':'bright_time_plates',
               'MWM_30min4':'bright_time_plates',
               'MWM3_sky':'bright_time_plates',
               'MWM4':'bright_time_plates',
               'OFFSET1':'eng_plates',
               'OFFSET2':'eng_plates',
               'RM':'dark_rm_plates',
               'RMv2':'dark_rm_plates',
               'RMv2-fewMWM':'dark_rm_plates'}
        try:
            designmode = p2d[programname]
        except:
            designmode = 'sdss5_plates'


    else:
        global chunkdata
        if chunkdata is None:
            chunkfile = ptt.join(os.getenv('PLATELIST_DIR'), 'platePlans.par')
            try:
                chunkdata = Table(yanny(chunkfile)['PLATEPLANS'])
            except:
                splog.info('Empty or missing platePlans.par file')
                chunkdata = None
        
        if chunkdata is not None:
            cinfo = chunkdata[np.where(chunkdata['plateid'] == int(plateid))[0]][0]
            programname = cinfo['programname']
            designmode = 'dark_'+programname+'_plates'
        else:
            designmode = 'dark_??_plates'

    #----------
    if not SOS:
        # Read calibObj or photoPlate photometry data
        programname = fibermap.meta['programname']
        ra_plate    = float(fibermap.meta['raCen'])
        dec_plate   = float(fibermap.meta['decCen'])
        fibermap    = calibrobj(fibermap, plateid, ra_plate, dec_plate, plates=plates,
                                legacy=legacy, no_db=no_db, mjd=mjd, indir=indir,
                                fast=fast, designmode=designmode, V_TARG=V_TARG,
                                RS_plan=fibermap.meta['platedesignversion'],
                                release=release, no_remote=no_remote)
        
    for col in ['org_fiberid','org_catid', 'carton', 'program_db',
                'll', 'bb', 'rr', 'stdflag']:
        if col in fibermap.colnames:
            fibermap.remove_column(col)
            
    fibermap = get_survey(fibermap, plates=plates, legacy=legacy)

    if legacy:
        fibermap['CatVersion']   = ''
    elif plates:
        fibermap['CatVersion']   = '0.0'
    fibermap.sort(['fiberId'])

    return(fibermap)


def get_reddening(search_table):
    global bayestar, no_bay
    global sfd, no_sfd
    global E3D2023, no_e3d
    global sd23, no_sd23
    
    splog.info('Calculating Reddening')
    search_table['EBV_rjce']=0.918 * ((search_table['TWOMASS_MAG'].data[:,1]) - (search_table['WISE_MAG'].data[:,1]) - 0.05)/0.302

    with HiddenPrints():
        gcord3d = SkyCoord(np.ma.filled(search_table['ll'].data)*u.deg,
                           np.ma.filled(search_table['bb'].data)*u.deg,
                           distance=np.ma.filled(search_table['rr'].data)*u.pc,
                           frame='galactic')

        gcord2d = SkyCoord(np.ma.filled(search_table['ll'].data)*u.deg,
                           np.ma.filled(search_table['bb'].data)*u.deg,
                           frame='galactic')
                        
        try:
            if 'bayestar' not in globals():
                bayestar = BayestarQuery(version='bayestar2015')
                no_bay=False
        except FileNotFoundError:
            try:
                import dustmaps.bayestar
                dustmaps.bayestar.fetch(version='bayestar2015')
                bayestar = BayestarQuery(version='bayestar2015')
                no_bay=False
            except ImportError: 
                no_bay=True
        try:
            if 'sd23' not in globals():
                sd23 = simple_dust_2023()
                no_sd23 = False
            no_e3d = True
        except:
            no_sd23 = True
            
            splog.warning('Simple_dust_2023 is unavailable... defaulting to Edenhofer2023')
            E3D2023_pars =    {'integrated':True, 'flavor': 'main'}
            E3D2023_2k_pars = {'integrated':True, 'flavor': 'less_data_but_2kpc'}
            try:
                if 'E3D2023' not in globals():
                    E3D2023 = E3D2023Query(**E3D2023_pars)
                    E3D2023 = 0
                    E3D2023_2k = E3D2023Query(**E3D2023_2k_pars)
                    E3D2023_2k = 0
                    no_e3d=False
            except FileNotFoundError:
                try:
                    import dustmaps.edenhofer2023
                    dustmaps.edenhofer2023.fetch(fetch_2kpc=True)
                    E3D2023 = E3D2023Query(**E3D2023_pars)
                    E3D2023 = 0
                    E3D2023_2k = E3D2023Query(**E3D2023_2k_pars)
                    E3D2023_2k = 0
                    no_e3d=False
                except ImportError:
                    no_e3d=True
        #no_e3d=True
        try:
            if 'sfd' not in globals():
                sfd = SFDQuery()
                no_sfd = False
        except FileNotFoundError:
            try:
                import dustmaps.sfd
                dustmaps.sfd.fetch()
                sfd = SFDQuery()
                no_sfd = False
            except ImportError:
                no_sfd=True
        #################################################################


        if not no_bay:
            splog.info('Getting Bayestar2015')
            ebv_bay = bayestar(gcord3d, mode='median')
            search_table['EBV_BAYESTAR15'] = ebv_bay
            
        if not no_sfd:
            splog.info('Getting SFD')
            search_table['SFD_EBV'] = sfd(gcord2d)

        if not no_sd23:
            splog.info('Getting SimpleDust2023')
            ebv_sd23 = sd23.query(gcord3d)
            search_table['EBV_SIMPLEDUST2023'] = ebv_sd23
            EBV_3D = ebv_sd23
            EBV_3DSRC = np.full(len(EBV_3D),'SimpleDust2023', dtype=object)
            if not no_bay:
                EBV_3D[np.where(np.isnan(EBV_3D))[0]] = ebv_bay[np.where(np.isnan(EBV_3D))[0]]
                EBV_3DSRC[np.where(np.isnan(EBV_3D))[0]] = 'bayestar15'
            search_table['EBV_3D'] = EBV_3D
            search_table['EBV_3DSRC'] = EBV_3DSRC

        elif not no_e3d:
            splog.info('Getting Edenhofer2023')
            E3D2023 = E3D2023Query(**E3D2023_pars)
            ebv_E3D2023 = E3D2023(gcord3d)
            E3D2023 = 0
            E3D2023_2k = E3D2023Query(**E3D2023_2k_pars)
            ebv_E3D2023_2k = E3D2023_2k(gcord3d)
            E3D2023_2k = 0
            lt2k = np.where(np.ma.filled(search_table['rr'].data)*u.pc > 1.25*u.kpc)[0]
            ebv_E3D2023[lt2k] = ebv_E3D2023_2k[lt2k]
            search_table['EBV_EDENHOFER2023'] = ebv_E3D2023
            
            EBV_3D = ebv_E3D2023
            EBV_3DSRC = np.full(len(EBV_3D),'edenhofer2023', dtype=object)
            EBV_3DSRC[lt2k] = 'edenhofer2023_2kpc'
            if not no_bay:
                EBV_3D[np.where(np.isnan(EBV_3D))[0]] = ebv_bay[np.where(np.isnan(EBV_3D))[0]]
                EBV_3DSRC[np.where(np.isnan(EBV_3D))[0]] = 'bayestar15'
            search_table['EBV_3D'] = EBV_3D
            search_table['EBV_3DSRC'] = EBV_3DSRC
            
        elif not no_bay:
            EBV_3D  = ebv_bay
            EBV_3DSRC = np.full(len(EBV_3D),'bayestar15', dtype=object)
            search_table['EBV_3D'] = EBV_3D
            search_table['EBV_3DSRC'] = EBV_3DSRC
    return(search_table)
    


def flag_too(search_table):
    if 'too' not in search_table.columns:
        return(search_table)
    if 'too_program' not in search_table.columns:
        program = search_table['program']
        toos = np.where((search_table['too'].data == 1) & (program == ''))[0]
        program[toos] = 'TOO'
        search_table['program'] = program
        return search_table

    program = search_table['program']

    too_program = search_table['too_program'].astype(str)  # Ensure all are strings
    too_program = np.char.upper(too_program).astype(object)  # Convert to uppercase and ensure it's an object array
    
    toos =  np.where((search_table['too'].data == 1))[0]

    toos = np.where((search_table['too'].data == 1) &
                    (too_program == ''))[0]
    too_program[toos] = 'TOO'
    program[toos] = too_program[toos]

    toos = np.where((search_table['too'] == 1) &
                    (~wwhere(too_program, '*TOO*')))[0]
    too_program[toos] = 'TOO_' + too_program[toos]
    program[toos] = too_program[toos]

    search_table['program'] = program
    search_table['too_program'] = too_program
    return(search_table)

def get_supplements(search_table, designID=None, rs_plan = None, fps=False, fast=False,
                    release='sdsswork', no_remote=False, db = True,
                    designmode=None,V_TARG='*'):

    with warnings.catch_warnings():
        warnings.simplefilter("error")
    
        dtypes = []
        for col in search_table.colnames:
            dtypes.append((col, search_table[col].dtype))
                          
        dtypes.extend([('CatVersion', object), ('carton', object), ('fieldCadence', object),
                       ('design_vers',object), ('DESIGN_MODE',object),
                       ('mapper', object), ('program_db', object), ('gaia_id', int), ('v05_rev_mag', bool),
                       ('EBV_rjce', float),('SFD_EBV',float), ('EBV_BAYESTAR15', float),
                       ('EBV_SIMPLEDUST2023',float),#('EBV_EDENHOFER2023',float),
                       ('EBV_3D',float), ('EBV_3DSRC', object),
                       ('ll', float), ('bb', float), ('rr', float),('SDSS_ID',int),
                       ('CATALOGID_V0',int),('CATALOGID_V0P5',int),('CATALOGID_V1',int)])
        if not fps:
            dtypes.extend([('parallax',float),('pmra',float),('pmdec',float)])
        data = Table(dtype=dtypes)
        for col in search_table.colnames:
            if len(search_table[col].shape) > 1:
                data[col] = Column(name = col, dtype=search_table[col].dtype, shape = (search_table[col].shape[1],))

        data['mag'] = Column(name='mag', dtype=float, shape=(5,))
        if fps is True:
            if designID is not None and rs_plan is not None:
                fieldCadence, ObsMode, designmode = get_FieldCadence(designID, rs_plan, db=db, release=release, 
                                                                     V_TARG=V_TARG, no_remote=no_remote)
                search_table['design_vers'] = rs_plan

            else:
                splog.info("No designID or rs_plan")
                fieldCadence = ''
                ObsMode = ''
                designmode = ''
                search_table['design_vers'] = ''

        else:
            fieldCadence = 'Plates'
            ObsMode = 'Plates'
            search_table['design_vers'] = 'plates_'+rs_plan
            if designmode is None:
                designmode = '??_plates'

        search_table['fieldCadence'] = fieldCadence
        search_table.meta['OBSMODE'] = ObsMode
        search_table.meta['DESIGN_MODE'] = designmode
        search_table['DESIGN_MODE'] = designmode

        for col in search_table.colnames:
            if search_table[col].dtype == float:
                search_table[col].fill_value = np.nan
            elif search_table[col].dtype == int:
                search_table[col].fill_value = -999
            elif search_table[col].dtype == object:
                search_table[col].fill_value = ''
            elif search_table[col].dtype == bool:
                search_table[col].fill_value = 0

        if fast is False:
            search_table = get_SDSSID(search_table, db=db, V_TARG=V_TARG, 
                                                 release=release, no_remote=no_remote)
            search_table, data = get_targetflags(search_table, data, db=db, V_TARG=V_TARG, 
                                                 release=release, no_remote=no_remote)
            search_table = get_AltCatids(search_table, db=db, V_TARG=V_TARG, 
                                                 release=release, no_remote=no_remote)


        search_table = get_mags_astrom(search_table, db = db, fps=fps, fast=fast,
                                       release=release, no_remote=no_remote,V_TARG=V_TARG)
        if (fps is True):
            if fast is False:
                search_table = get_CartonInfo(search_table, db=db, V_TARG=V_TARG, 
                                                 release=release, no_remote=no_remote)
            search_table = target_tab_correction(search_table, db=db, V_TARG=V_TARG, 
                                                 release=release, no_remote=no_remote)
        calc_dist=True
        if calc_dist is True:
            gcord = SkyCoord(search_table['ra'].data*u.deg, search_table['dec'].data*u.deg).transform_to('galactic')
            search_table['ll']=gcord.l.deg
            search_table['bb']=gcord.b.deg
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                search_table['rr']=Distance(parallax=search_table['parallax'].data*u.mas,allow_negative=True).value
        search_table = get_reddening(search_table)
    
        search_table = flag_too(search_table)


        for col in ['parallax','pmra','pmdec','optical_prov',
                    'mag_u','mag_g','mag_r','mag_i','mag_z',
                    'mag_h','gaia_g','gaia_bp','gaia_rp']:
            if col not in search_table.colnames:
                continue
            cval = []
            for i, row in enumerate(search_table):
                if col not in search_table.colnames:
                    continue
                if ((search_table[col].dtype == float) or (search_table[col].dtype == int)):
                    if np.ma.is_masked(search_table[col].data[i]):
                        indx = np.where(search_table['fiberId'] == search_table[i]['fiberId'])[0]
                        if len(indx) > 0:
                            val = np.nan #search_table[col][indx].data
                        cval.append(val)
                    elif np.isnan(search_table[col].data[i]):
                        indx = np.where(search_table['fiberId'] == search_table[i]['fiberId'])[0]
                        if len(indx) > 0:
                            val =  np.nan #search_table[col][indx].data
                        cval.append(val)
                    else:
                        cval.append(search_table[col].data[i])
                else:
                    if np.ma.is_masked(search_table[col].data[i]):
                        indx = np.where(search_table['fiberId'] == search_table[i]['fiberId'])[0]
                        if len(indx) > 0:
                            val = ''#search_table[col][indx].data
                        cval.append(val)
                    elif np.ma.getdata(search_table[col].data, subok=False)[i] == '':
                        indx = np.where(search_table['fiberId'] == search_table[i]['fiberId'])[0]
                        if len(indx) > 0:
                            val = ''#search_table[col][indx].data
                        cval.append(val)
                    else:
                        cval.append(search_table[col].data[i])
            cval = np.asarray(cval)
            if search_table[col].dtype == object:
                cval = cval.astype(object)
            for i, val in enumerate(cval):
                if cval[i] == -999.:
                    cval[i] = np.nan



            if search_table[col].dtype == object:
                cval=cval.astype(object)
                search_table.remove_column(col)
                search_table.add_column(cval, name=col)
            else:
                search_table.remove_column(col)
                search_table.add_column(cval, name=col)

        if fps:
            for i, row in enumerate(data):
                if (np.ma.is_masked(data['mapper'].data[i])) or (np.ma.getdata(data['mapper'].data, subok=False)[i] == ''):
                    if row['program'] == 'open_fiber':
                        search_table[i]['mapper'] = 'open_fiber'
                    elif 'ops' in row['program']:
                        search_table[i]['mapper'] = 'ops'
                    else:
                        search_table[i]['mapper'] = ''


        search_table['catalogid'] = search_table['catalogid'].astype(object)

        for col in search_table.colnames:
            if not ((search_table[col].dtype  in [float, int, np.dtype('int16'), np.dtype('float32'), bool])):
                search_table[col] = search_table[col].astype(object)
        for col in data.colnames:
            if not ((data[col].dtype  in [float, int, np.dtype('int16'), np.dtype('float32'), bool])):
                data[col] = data[col].astype(object)

        search_table = vstack([data,search_table])
        search_table = search_table[data.colnames]

        ops_std = np.zeros(len(search_table), dtype=bool)
        ops_std[np.where((search_table['program'].data == b'ops_std') | (search_table['program'].data == 'ops_std'))[0]] = 1
        search_table['stdflag'] = ops_std

        for col in search_table.colnames:
            if search_table[col].dtype == float:
                search_table[col].fill_value = np.nan
            elif search_table[col].dtype == int:
                search_table[col].fill_value = -999
            elif search_table[col].dtype == object:
                search_table[col].fill_value = ''
            elif search_table[col].dtype == bool:
                search_table[col].fill_value = 0

        search_table=search_table.filled()

    return(search_table)





    
