#!/usr/bin/env python3


try:
    from load_module import load_module, load_env
except:
    from load_module import load_module
    def load_env(key, default=None):
        val = getenv(key)
        if val is None:
            val = default
            if val is None:
                print('ERROR: '+key+' is not set')
                exit()
        return(val)
import argparse
from sdssdb.peewee.sdss5db.targetdb import database

test = database.set_profile('pipelines')

from sdssdb.peewee.sdss5db.targetdb import CartonToTarget, Carton, Version, Mapper, Target
from astropy.table import Table, Column
from astropy.io import fits
from os import getenv, remove
import os.path as ptt
from glob import glob
import numpy as np
from tqdm.auto import tqdm
import time
from readfibermaps import get_targetflags
from splog import Splog
import shutil

splog = Splog()
def get_Targeting_file(run2d, boss_spectro_redux=getenv('BOSS_SPECTRO_REDUX')):
    return(ptt.join(boss_spectro_redux, run2d, f'spTargeting-{run2d}.fits.gz'))


def update_Targeting_flags(run2d, boss_spectro_redux, schema = None, clobber=False,
                           nobackup=False):
    spall_file = ptt.join(boss_spectro_redux, run2d, f'spAll-lite-{run2d}.fits.gz')
    Targeting_file = get_Targeting_file(run2d, boss_spectro_redux=boss_spectro_redux)
    
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

    run2d_dir = ptt.join(boss_spectro_redux, run2d)
    files = ['spAll-lite-{run2d}{schema}.fits.gz','spAll-{run2d}{schema}.fits.gz']
    if schema is not None:
        custom = schema
        schema = [None]
        schema.extend(custom)
    else:
        schema = [None]
    for epoch in ['','epoch']:
        for sc in schema:
            sc = '' if sc is None else '-'+sc
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
                        idx = np.searchsorted(search_table['SDSS_ID'],hdul[1].data['SDSS_ID'])
                        nf  = np.where(hdul[1].data['SDSS_ID'] != -999)[0]
    
                        spall_shape = hdul[1].data['SDSS5_TARGET_FLAGS'].shape[1]
                        update_shape = search_table['SDSS5_TARGET_FLAGS'].shape[1]
                        if spall_shape < update_shape:
                            splog.info(f'Padding SDSS5_TARGET_FLAGS in {ptt.basename(ff)}')
                            coldat = hdul[1].data['SDSS5_TARGET_FLAGS']
                            pad = update_shape - spall_shape
#                            hdul[1].data['SDSS5_TARGET_FLAGS'] = np.pad(coldat,[(0,0),(pad,0)], mode = 'constant', constant_values= 0)
#                            hdul[1].data['SDSS5_TARGET_FLAGS'] = fits.Column(name='SDSS5_TARGET_FLAGS', format=f'{update_shape}B',
#                                                                            array=np.pad(coldat,[(0,0),(pad,0)], mode = 'constant', constant_values= 0))
                            new_col = fits.Column(name='SDSS5_TARGET_FLAGS', format=f'{update_shape}B',
                                                  array=np.pad(coldat,[(0,0),(pad,0)], mode = 'constant', constant_values= 0))
                            new_col = fits.ColDefs([new_col])
                            old_cols = hdul[1].columns
                            colidx = np.where(np.asarray(old_cols.names) == 'SDSS5_TARGET_FLAGS')[0][0]
                            
                            hdul[1] = fits.BinTableHDU.from_columns(old_cols[0:colidx-1]+new_col+old_cols[colidx+1:], 
                                                                    header=hdul[1].header)
                        hdul[1].data['SDSS5_TARGET_FLAGS'][nf] = search_table['SDSS5_TARGET_FLAGS'][idx][nf] 
                        #hdul[1].data['SDSSC2BV'] = search_table['SDSSC2BV'][idx][nf]
                        hdul[0].header['SDSSC2BV'] = search_table['SDSSC2BV'][0]

    
    

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot QA')
    parser.add_argument('-m','--mod',   required=False,       help='bhm module to run', default = 'bhm/master',)
    parser.add_argument('--clobber',    action='store_true',  help='Clobber spTargeting file')
    parser.add_argument('--custom',     nargs='*',            help='List of name of custom coadd schema')
    parser.add_argument('--nobackup',   action='store_true',  help='Skip backup of existing')
    args = parser.parse_args()

    module = load_module()
    module('purge')
    module('load', args.mod)
    run2d = load_env('RUN2D')
    boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    
    update_Targeting_flags(run2d, boss_spectro_redux, schema=args.custom, clobber=args.clobber,
                           nobackup=args.nobackup)
