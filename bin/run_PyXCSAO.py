#!/usr/bin/env python

import os
import numpy as np
try:
    import pyxcsao
    from pyxcsao.crosscorrelate import PyXCSAO
except:
    print('WARNING: pyxcsao is not installed')
from astropy.table import Table
import pandas as pd
import argparse
import sys
from astropy.io import fits
from astropy.table import Table
from platform import python_version
from datetime import datetime
"""
run_PyXCSAO.py:

run_PyXCSAO.py runs PyXCSAO for a full spField-*****-*****.fits file using the phoenix_full1 template grid
The input files can be either normal or gz files.
"""

class LogFile:
    def __init__(self,logfile):
        self.splog = open(logfile, 'w')
    def write(self,logline):
        logline = 'run_PyXCSAO: '+ logline
        logline=logline.replace('\n','')
        logline=logline.replace('\r','')
        self.splog.write(logline+'\n')
        print(logline)
    def close(self):
        self.splog.write('\n')
        self.splog.close()
        
def load_boss_field(name):
    with fits.open(name) as hdul:
        flux = hdul[0].data
        PlugMap = hdul['PLUGMAP'].data
        hdr = hdul[0].header
        la = 10**(np.arange(hdr['NAXIS1'])*hdr['CD1_1']+hdr['CRVAL1'])
    return la, flux, PlugMap, hdr

def get_fiber(flux, PlugMap, hdr, i):
    meta={}
    meta['ra']=PlugMap['RA'][i]
    meta['dec']=PlugMap['DEC'][i]
    meta['objid']=PlugMap['CATALOGID'][i]
    #meta['objid']=PlugMap['OBJID'][i]
    meta['program']=PlugMap['PROGRAM'][i]
    meta['objtype']=PlugMap['OBJTYPE'][i]
    meta['SOURCETYPE']=PlugMap['SOURCETYPE'][i]
    
    if 'PLATEID' in hdr: meta['plate']=hdr['PLATEID']
    elif 'FIELDID' in hdr: meta['plate']=hdr['FIELDID']
    else: meta['plate']=np.nan
    
    meta['mjd']=hdr['MJD']
    meta['fiber']=PlugMap['FIBERID'][i]
    #meta['snr']=PlugMap['SN_MEDIAN_ALL'][i]
    meta['snr']=np.nan
    meta['firstcarton']=PlugMap['FIRSTCARTON'][i]
    meta['parallax']=PlugMap['GAIA_PARALLAX'][i]
    meta['pmra']=PlugMap['GAIA_PMRA'][i]
    meta['pmdec']=PlugMap['GAIA_PMDEC'][i]
    
    meta['SFD_EBV']=PlugMap['SFD_EBV'][i]

    meta['sdss_u']=PlugMap['MAG'][i][0]
    meta['sdss_g']=PlugMap['MAG'][i][1]
    meta['sdss_r']=PlugMap['MAG'][i][2]
    meta['sdss_i']=PlugMap['MAG'][i][3]
    meta['sdss_z']=PlugMap['MAG'][i][4]
    meta['gaia_G']=PlugMap['GAIA_G'][i]
    meta['BP']=PlugMap['GAIA_BP'][i]
    meta['RP']=PlugMap['GAIA_RP'][i]
    meta['J']=PlugMap['TWOMASS_MAG'][i][0]
    meta['H']=PlugMap['TWOMASS_MAG'][i][1]
    meta['K']=PlugMap['TWOMASS_MAG'][i][2]
    
    return flux[i,:], meta

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Runs eBOSS RVs')

    parser.add_argument('fitsfile', type=str, help='fits file')
    parser.add_argument('--run1d', '-r', type=str, help='run1d name',
                        default=os.getenv('RUN1D'))
    args = parser.parse_args()

    platemjd=os.path.basename(args.fitsfile)[8:19]


    # Open a file with access mode 'a'
    logfile = args.run1d+'/spXCSAO-' + platemjd + '.log'
    splog = LogFile(logfile)
    t0=datetime.now()
    splog.write('Log file ' +logfile +' opened '+t0.strftime("%a %b %d %H:%M:%S %Y"))
    splog.write('UNAME: '+ os.popen('uname -a').read())
    try: splog.write('DISPLAY= ' + os.getenv('DISPLAY'))
    except: splog.write('DISPLAY= ')
    splog.write('idlspec2d version ' + os.popen('idlspec2d_version').read())
    splog.write('idlutils version ' + os.popen('idlutils_version').read())
    splog.write('Python version ' + python_version())
    try:
        #splog.write('pyxcsao version ' + pyxcsao.__path__[0].split('/')[-2][:-4])
        splog.write('pyxcsao version ' +os.getenv('PYXCSAO_VER'))
        c=PyXCSAO(st_lambda=5000,end_lambda=10000)
        templates=os.getenv('PYXCSAO_DIR')+'/../grids/phoenix_full1.p'
        c.add_grid(grid_pickle=templates)

        best=[]
        la, flux, PlugMap, hdr = load_boss_field(args.fitsfile)
        for i in range(hdr['NAXIS2']):
            splog.write('Calculating rv for fiber '+str(i))
            flux_one, meta = get_fiber(flux, PlugMap, hdr, i)
            try:
                c.add_spectrum(flux_one,i=i,meta=meta,laname=la, data_class='user')
                out = c.run_XCSAO_optimized().copy()
                best.append(out)
            except:
                best.append(meta)
        df=pd.DataFrame(best)

        df['objid']=df['objid'].astype(str)
        test=Table.from_pandas(df)
    
        test.write(args.run1d+'/spXCSAO-' + platemjd + '.fits',overwrite=True)

        splog.write('CPU time to compute RVs = '+ str(datetime.now()-t0))
    except: splog.write('WARNING: Python package pyxcsao not installed\n')
    splog.close()
