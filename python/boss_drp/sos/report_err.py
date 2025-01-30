#!/usr/bin/env python3

from boss_drp.sos.run_log2html import run_soslog2html
from boss_drp.utils.lock import unlock, lock

import os
import os.path as ptt
from astropy.io import fits
from astropy.table import Table, unique
import traceback
import time


def report(FitsName, cams, obs, mjd, message, designMode = 'unknown'):
    i = 0
    if obs.lower() == 'apo':
        sos_data = 'BOSS_SPECTRO_DATA_N'
    else:
        sos_data = 'BOSS_SPECTRO_DATA_S'
    hdr = fits.getheader(ptt.join(os.getenv(sos_data),f'{mjd}',FitsName))
    logfile = ptt.join(os.getenv('BOSS_SPECTRO_REDUX'), f'logfile-{mjd}.fits')
    
    if lock(logfile,niter=10, pause=2):
        try:
            with fits.open(logfile) as hdul:
                if isinstance(hdul[5], fits.BinTableHDU):
                    try:
                        log = Table.read(logfile,hdu=5)
                        for col in ['FILENAME','CARTID','CAMERA','EXPNUM','TEXT','DESIGNMODE']:
                            print(col)
                            log[col] = log[col].astype(object)
                    except Exception as e:
                        print(e)
                        log = Table(names=['FILENAME','MJD','CONFIG','FIELD','CARTID','DESIGNMODE','EXPNUM','CAMERA','TEXT'],
                                    dtype=[object,int,int,int,object,object,object,object,object])
                else:
                    log = Table(names=['FILENAME','MJD','CONFIG','FIELD','CARTID','DESIGNMODE','EXPNUM','CAMERA','TEXT'],
                               dtype=[object,int,int,int,object,object,object,object,object])
            try:
                log = Table.read(logfile,hdu=5)
                for col in ['FILENAME','CARTID','CAMERA','EXPNUM','TEXT']:
                    log[col] = log[col].astype(object)
            except Exception as e:
                print(e)
                log = Table(names=['FILENAME','MJD','CONFIG','FIELD','CARTID','EXPNUM','CAMERA','TEXT'],
                            dtype=[object,int,int,int,object,object,object,object])

            log.add_row([hdr['FILENAME'].encode('utf-8'),
                         hdr['MJD'],
                         hdr['CONFID'],
                         hdr['FIELDID'],
                         hdr['CARTID'].strip().encode('utf-8'),
                         designMode.strip().encode('utf-8'),
                         str(hdr['EXPOSURE']).zfill(8).strip().encode('utf-8'),
                         cams.strip().encode('utf-8'),
                         message.encode('utf-8')])
            log = unique(log, keys=['FILENAME','CAMERA','TEXT'], keep='last')

            for col in ['FILENAME','CARTID','CAMERA','EXPNUM','TEXT','DESIGNMODE']:
                log[col] = log[col].astype(str)
            if ptt.exists(logfile):
                with fits.open(logfile, mode='update', output_verify="silentfix") as hdul:
                    hdul[5] = fits.table_to_hdu(log)
                    hdul.flush()
            else:
                hdul = fits.HDUList()
                hdul.append(fits.PrimaryHDU()) # 0
                hdul.append(fits.ImageHDU()) # bias/dark
                hdul.append(fits.ImageHDU()) # flat
                hdul.append(fits.ImageHDU()) # arc
                hdul.append(fits.ImageHDU()) # science
                hdul.append(fits.table_to_hdu(log)) # Messages
                hdul.writeto(logfile)
        except Exception as e:
            tb_str = traceback.format_exception(type(e), e, e.__traceback__)
            print("".join(tb_str))
            e = tb_str = None
        finally:
            unlock(logfile)
    else:
        print("Could not acquire lock. Exiting.")

    run_soslog2html(logfile, mjd, obs)
    return
