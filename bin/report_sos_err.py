import os
import os.path as ptt
from astropy.io import fits
from astropy.table import Table, unique
import time
import traceback
from run_soslog2html import run_soslog2html

def report(FitsName, cams, obs, mjd, message):
    i = 0
    if obs.lower() == 'apo':
        sos_data = 'BOSS_SPECTRO_DATA_N'
    else:
        sos_data = 'BOSS_SPECTRO_DATA_S'
    hdr = fits.getheader(ptt.join(os.getenv(sos_data),f'{mjd}',FitsName))
    logfile = ptt.join(os.getenv('BOSS_SPECTRO_REDUX'), f'logfile-{mjd}.fits')
    while True:
        if i != 0: time.sleep(2)
        i = i+1
        if i > 10: break
        linked = False
        try:
            try:
                os.symlink(logfile, logfile+'.lock')
                linked = True
            except Exception as e:
                print(e)
                continue
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
                         str(hdr['EXPOSURE']).zfill(8).strip().encode('utf-8'),
                         cams.strip().encode('utf-8'),
                         message.encode('utf-8')])
            log = unique(log, keys=['FILENAME','CAMERA','TEXT'], keep='last')

            for col in ['FILENAME','CARTID','CAMERA','EXPNUM','TEXT']:
                log[col] = log[col].astype(str)
            if ptt.exists(logfile):
                with fits.open(logfile, mode='update') as hdul:
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
            os.unlink(logfile+'.lock')
            break
        except Exception as e:
            tb_str = traceback.format_exception(etype=type(e), value=e, tb=e.__traceback__)
            print("".join(tb_str))
            e = tb_str = None
            if linked:
                print('failed... unlocking')
                os.unlink(logfile+'.lock')
            continue
    
    run_soslog2html(logfile, mjd, obs)
    return