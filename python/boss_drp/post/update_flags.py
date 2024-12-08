#!/usr/bin/env python3
from boss_drp.prep.readfibermaps import get_targetflags
from boss_drp.utils import (load_env, Splog)

splog = Splog()

try:
    from sdssdb.peewee.sdss5db.targetdb import database
    test = database.set_profile(load_env('DATABASE_PROFILE', default='pipelines'))
    from sdssdb.peewee.sdss5db.targetdb import CartonToTarget, Carton, Version, Mapper, Target
    nodb=False
except:
    if load_env('DATABASE_PROFILE', default='pipelines').lower() in ['pipelines','operations']:
        splog.log('ERROR: No SDSSDB access')
        exit()
    else:
        splog.log('ERROR: No SDSSDB access')
    nodb= True

import argparse
from astropy.table import Table, Column, join
from astropy.io import fits
from os import getenv, remove
import os.path as ptt
from glob import glob
import numpy as np
import time
import shutil




def get_Targeting_file(run2d, boss_spectro_redux=getenv('BOSS_SPECTRO_REDUX')):
    return(ptt.join(boss_spectro_redux, run2d, 'summary','daily', f'spTargeting-{run2d}.fits.gz'))


def update_Targeting_flags(run2d, boss_spectro_redux, schema = None, clobber=False,
                           nobackup=False,nrows=None, spall_tab=None, idcol='CATALOGID',
                           legacy=False):
    if nodb:
        exit()
    spall_file = ptt.join(boss_spectro_redux, run2d, 'summary','daily', f'spAll-lite-{run2d}.fits.gz')
    Targeting_file = get_Targeting_file(run2d, boss_spectro_redux=boss_spectro_redux)
    if not ptt.exists(Targeting_file):
        clobber = True
    if clobber:
        splog.info('Determining Updated Flags')
        spall = fits.getdata(spall_file, 1)
        sdssids = np.unique(spall['SDSS_ID'].data)
        sdssids = np.setdiff1d(sdssids,np.array([-999,0]))
        search_table = Table([sdssids], names = ['SDSS_ID'])
        search_table, _junk= get_targetflags(search_table, None)
        for col in ['SDSS5_TARGET_FLAGS','SDSSC2BV']:
            search_table[col].fill_value = 0
            #search_table[col][search_table[col].data.astype(object) == ''] = 0
            
            try:
                data = search_table[col].data.filled()
            except:
                data = search_table[col].data
            search_table[col] = data.astype('uint8')
        search_table.write(get_Targeting_file(run2d, boss_spectro_redux), overwrite=True)
    else:
        if ptt.exists(Targeting_file):
            splog.info('Reading: '+Targeting_file)
            search_table = Table.read(Targeting_file)
        else:
            raise Exception(Targeting_file+' is missing')

    run2d_dir = ptt.join(boss_spectro_redux, run2d, 'summary')
    files = ['spAll-lite-{run2d}{schema}.fits.gz','spAll-{run2d}{schema}.fits.gz']
        
    epochs = ['daily','epoch']
    if schema is not None:
        epochs.extend(schema)
    for epoch in epochs:
        if schema is not None:
            sc = '' if epoch in ['daily'] else '-'+epoch
        else:
            sc = '' if epoch in ['daily'] else '-'+epoch
        for f in files:
            ff = ptt.join(run2d_dir,epoch,f.format(run2d=run2d, schema = sc))
            
            if ptt.exists(ff):
                if not nobackup:
                    if ptt.exists(ff+'.bak'):
                        splog.info(f'Removing {ptt.basename(ff)}.bak')
                        remove(ff+'.bak')
                    try:
                        splog.info(f'Saving a backup of {ff} to {ff}.bak')
                        shutil.copy(ff,ff+'.bak')
                    except OSError as exc:
                        raise OSError(f"Failed to save backup to destination {ff}") from exc
                with fits.open(ff, mode='update', lazy_load_hdus=True, save_backup=False) as hdul:
                    splog.info(f'Updating {ff}')
                    nf  = np.where(hdul[1].data['SDSS_ID'] != -999)[0]
                    idx = np.searchsorted(search_table['SDSS_ID'],hdul[1].data['SDSS_ID'][nf])

                    spall_shape = hdul[1].data['SDSS5_TARGET_FLAGS'].shape[1]
                    update_shape = search_table['SDSS5_TARGET_FLAGS'].shape[1]
                    if spall_shape < update_shape:
                        splog.info(f'Padding SDSS5_TARGET_FLAGS in {ptt.basename(ff)}')
                        coldat = hdul[1].data['SDSS5_TARGET_FLAGS']
                        pad = update_shape - spall_shape
                        new_col = fits.Column(name='SDSS5_TARGET_FLAGS', format=f'{update_shape}B',
                                              array=np.pad(coldat,[(0,0),(0,pad)], mode = 'constant', constant_values= 0))
                        new_col = fits.ColDefs([new_col])
                        old_cols = hdul[1].columns
                        colidx = np.where(np.asarray(old_cols.names) == 'SDSS5_TARGET_FLAGS')[0][0]
                        
                        hdul[1] = fits.BinTableHDU.from_columns(old_cols[0:colidx-1]+new_col+old_cols[colidx+1:],
                                                                header=hdul[1].header)
                    elif spall_shape > update_shape:
                        pad = spall_shape - update_shape
                        coldat = search_table['SDSS5_TARGET_FLAGS'].data

                        search_table['SDSS5_TARGET_FLAGS'] = np.pad(coldat, [(0,0),(0,pad)],
                                                                    mode = 'constant',
                                                                    constant_values= 0)


                    search_table_j = join(Table([hdul[1].data['SDSS_ID'], range(len(hdul[1].data))], names=['SDSS_ID','order']), search_table, keys='SDSS_ID',join_type='left')
                    search_table_j.sort('order')
                    hdul[1].data['SDSS5_TARGET_FLAGS'][nf] = search_table_j['SDSS5_TARGET_FLAGS'][nf]
                    hdul[0].header['SDSSC2BV'] = search_table_j['SDSSC2BV'][0]

