#!/usr/bin/env python

import os
import numpy as np
try:
    import pyxcsao
    from pyxcsao.crosscorrelate import PyXCSAO
except:
    print('WARNING: pyxcsao is not installed')
import pandas as pd
import argparse
import sys
from astropy.io import fits
from astropy.table import Table
from platform import python_version
from datetime import datetime
from splog import Splog

"""
run_PyXCSAO.py:

run_PyXCSAO.py runs PyXCSAO for a full spField-*****-*****.fits file using the phoenix_full1 template grid
The input files can be either normal or gz files.
"""
splog = Splog()
        
def load_boss_field(name):
    with fits.open(name) as hdul:
        flux = hdul[0].data
        PlugMap = hdul['PLUGMAP'].data
        hdr = hdul[0].header
        la = 10**(np.arange(hdr['NAXIS1'])*hdr['CD1_1']+hdr['CRVAL1'])
    return la, flux, PlugMap, hdr

def get_fiber(flux, PlugMap, hdr, i):
    meta={}
    if 'ra' in PlugMap.columns.names:
        meta['ra']=PlugMap['ra'][i]
    elif 'RA' in  PlugMap.columns.names:
        meta['ra']=PlugMap['RA'][i]
    else:
        meta['ra']=PlugMap['FIBER_RA'][i]
        
    if 'dec' in PlugMap.columns.names:
        meta['dec']=PlugMap['dec'][i]
    elif 'DEC' in PlugMap.columns.names:
        meta['dec']=PlugMap['DEC'][i]
    else:
        meta['dec']=PlugMap['FIBER_DEC'][i]
        
    if 'coord_epoch' in PlugMap.columns.names:
        meta['coord_epoch'] = PlugMap['coord_epoch'][i]
    elif 'COORD_EPOCH' in PlugMap.columns.names:
        meta['coord_epoch'] = PlugMap['COORD_EPOCH'][i]
    else:
        meta['coord_epoch'] = 2000


    meta['objid']=PlugMap['CATALOGID'][i]
    #meta['objid']=PlugMap['OBJID'][i]
    meta['program']=PlugMap['PROGRAM'][i]
    meta['objtype']=PlugMap['OBJTYPE'][i]
    if 'CATEGORY' in PlugMap.columns.names:
        meta['SOURCETYPE']=PlugMap['CATEGORY'][i]
    elif 'SOURCETYPE' in PlugMap.columns.names:
        meta['SOURCETYPE']=PlugMap['SOURCETYPE'][i]
    else:
        meta['SOURCETYPE']=PlugMap['OBJTYPE'][i]


    if 'FIELDID' in hdr:
        meta['FIELDID']=hdr['FIELDID']
    elif 'PLATEID' in hdr:
        meta['FIELDID']=hdr['PLATEID']
    elif 'FIELDID' in hdr:
        meta['FIELDID']=hdr['FIELDID']
    elif 'COADD' in hdr:
        meta['FIELDID']=hdr['COADD']
    else:
        meta['FIELDID']=np.nan
    if meta['FIELDID'].strip() == '':
        if 'PLATEID' in hdr:
            meta['FIELDID']=hdr['PLATEID']
    
    meta['mjd']=hdr['MJD']
    meta['TARGET_INDEX']=PlugMap['TARGET_INDEX'][i]
    if 'FIBERID_LIST' in PlugMap.names:
        meta['FIBERID_LIST']=PlugMap['FIBERID_LIST'][i]
    #meta['fiber']=PlugMap['FIBERID'][i]
    #meta['snr']=PlugMap['SN_MEDIAN_ALL'][i]
    meta['snr']=np.nan
    meta['firstcarton']=PlugMap['FIRSTCARTON'][i]
    meta['SDSS_ID']=PlugMap['SDSS_ID'][i]
    meta['CATALOGID']=PlugMap['ICATALOGID'][i]
    meta['firstcarton']=PlugMap['FIRSTCARTON'][i]

    meta['parallax']=PlugMap['PARALLAX'][i]
    meta['pmra']=PlugMap['PMRA'][i]
    meta['pmdec']=PlugMap['PMDEC'][i]
    
    meta['EBV']=PlugMap['EBV'][i]

    meta['sdss_u']=PlugMap['MAG'][i][0]
    meta['sdss_g']=PlugMap['MAG'][i][1]
    meta['sdss_r']=PlugMap['MAG'][i][2]
    meta['sdss_i']=PlugMap['MAG'][i][3]
    meta['sdss_z']=PlugMap['MAG'][i][4]
    if 'GAIA_G' in PlugMap.names:
        meta['gaia_G']=PlugMap['GAIA_G'][i]
    elif 'GAIA_G_MAG' in PlugMap.names:
        meta['gaia_G']=PlugMap['GAIA_G_MAG'][i]
    else:
        meta['gaia_G']=np.NaN
        
    if 'GAIA_BP' in PlugMap.names:
        meta['BP']=PlugMap['GAIA_BP'][i]
    elif 'BP_MAG' in PlugMap.names:
        meta['BP']=PlugMap['BP_MAG'][i]
    else:
        meta['BP']=np.NaN
        
    if 'GAIA_RP' in PlugMap.names:
        meta['RP']=PlugMap['GAIA_RP'][i]
    elif 'RP_MAG' in PlugMap.names:
        meta['RP']=PlugMap['RP_MAG'][i]
    else:
        meta['RP']=np.NaN
    meta['J']=PlugMap['TWOMASS_MAG'][i][0]
    meta['H']=PlugMap['TWOMASS_MAG'][i][1]
    meta['K']=PlugMap['TWOMASS_MAG'][i][2]
    
    return flux[i,:], meta


def run_PyXCSAO(fitsfile,run1d = os.getenv('RUN1D'), epoch=False):

    platemjd='-'.join(os.path.basename(fitsfile).split('-')[1:]).split('.')[0]

    # Open a file with access mode 'a'
    logfile = run1d+'/spXCSAO-' + platemjd + '.log'
    splog.open(logfile=logfile)
    t0=datetime.now()
    splog.log('Log file ' +logfile +' opened '+t0.strftime("%a %b %d %H:%M:%S %Y"))
    splog.log('UNAME: '+ os.popen('uname -a').read())
    splog.log('DISPLAY= ' + os.getenv('DISPLAY', default=''))
    splog.log('idlspec2d version ' + os.popen('idlspec2d_version').read())
    splog.log('idlutils version ' + os.popen('idlutils_version').read())
    splog.log('Python version ' + python_version())
    if True:
    #try:
        splog.log('pyxcsao version ' +os.getenv('PYXCSAO_VER'))
        c=PyXCSAO(st_lambda=5000,end_lambda=10000)
        templates = os.path.join(os.getenv('PYXCSAO_DIR'), '..','grids','phoenix_full1.p')
        c.add_grid(grid_pickle=templates)

        best=[]
        la, flux, PlugMap, hdr = load_boss_field(fitsfile)
        for i in range(hdr['NAXIS2']):
            splog.log('Calculating rv for TARGET_INDEX '+str(i+1)+'/'+str(hdr['NAXIS2']))
            flux_one, meta = get_fiber(flux, PlugMap, hdr, i)
            try:
                c.add_spectrum(flux_one,i=i,meta=meta,laname=la, data_class='user')
                out = c.run_XCSAO_optimized(optimized_for_boss=True).copy()
                best.append(out)
            except:
                best.append(meta)
        df=pd.DataFrame(best)

        df['objid']=df['objid'].astype(str)
        df.drop(columns=['snr', 'plate', 'fiber'], inplace=True, errors='ignore')
        test=Table.from_pandas(df)
    
        test.write(run1d+'/spXCSAO-' + platemjd + '.fits',overwrite=True)

        splog.log('CPU time to compute RVs = '+ str(datetime.now()-t0))
    #except: splog.log('WARNING: Failed run of pyXCSAO\n')
    splog.close()




if __name__ == '__main__' :

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Runs eBOSS RVs')

    parser.add_argument('fitsfile', type=str, help='fits file')
    parser.add_argument('--run1d', '-r', type=str, help='run1d name',
                        default=os.getenv('RUN1D'))
    parser.add_argument('--epoch', help='run for epoch Coadds', action='store_true')
    args = parser.parse_args()

    run_PyXCSAO(args.fitsfile, run1d = args.run1d, epoch = False)
