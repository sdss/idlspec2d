from boss_drp import idlspec2d_dir

from astropy.io import fits
from pydl.pydlutils.yanny import read_table_yanny
import numpy as np
import os
import warnings

def refine_opfiber(spFlat, precision = 3):
    cart = fits.getheader(spFlat)['CARTID']
    cam  = fits.getheader(spFlat)['CAMERAS']
    print('MJD: ',fits.getheader(spFlat)['MJD'])
    
    opFibers = read_table_yanny(os.path.join(idlspec2d_dir, 'opfiles', 'opFibersFPS.par'), 'FIBERPARAM')
    opFibers = opFibers[(opFibers['cartid'] == cart) & (opFibers['camname'] ==cam)]
    opFibers= opFibers[opFibers['mjd'] == opFibers['mjd'].max()]

    xcen = fits.getdata(spFlat,5)
    xcen = np.mean(xcen[:,xcen.shape[1]//2-5:xcen.shape[1]//2+5],axis=1)

    bundlegaps = []
    fiberspacing = []
    sfiber = 0
    for i, nfib in enumerate(opFibers['bundlefibers'][0]):
        if i == 0:
            bundlegaps.append(xcen[0])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            fspace = np.mean(np.diff(xcen[sfiber:sfiber+nfib]))
            if not np.isfinite(fspace):
                if i == 0:
                    fspace = 0
                elif nfib == 1:
                    fspace = np.mean(fiberspacing)
                else:
                    fspace = 0
            if i < opFibers['nbundles'][0] - 1:
                try:
                    bgap = [np.diff(xcen[sfiber+nfib-1: sfiber+nfib+1])[0] - fspace]
                except:
                    bgap = [np.NaN]
            elif i == opFibers['nbundles'][0] - 1:
                bgap = [0.0]
            else:
                bgap = []
        bundlegaps.extend(bgap)
        fiberspacing.append(fspace)
        sfiber = sfiber+nfib

    format_str = f".{precision}f"
    fiberspacing = '{ '+" ".join(f"{x:{format_str}}" for x in fiberspacing)+' }'
    bundlegaps = '{ '+" ".join(f"{x:{format_str}}" for x in bundlegaps)+' }'
    
    bundlefibers = '{ '+" ".join(f"{x:d}" for x in opFibers['bundlefibers'][0])+' }'
    deadfibermask = '{ '+" ".join(f"{x:d}" for x in opFibers['deadfibermask'][0])+' }'

    print('--------------------- New version ----------------------')

    print(f"FIBERPARAM {cart} {cam} {opFibers['mjd'][0]} 70.0000 {np.sum(opFibers['bundlefibers'][0])} "+
          f"{opFibers['nbundles'][0]} {fiberspacing} {bundlegaps} {bundlefibers} {deadfibermask}")
 
    print('------------------ Previous version --------------------')
    fiberspacing = '{ '+" ".join(f"{x:{format_str}}" for x in opFibers['fiberspace'][0])+' }'
    bundlegaps = '{ '+" ".join(f"{x:{format_str}}" for x in opFibers['bundlegap'][0])+' }'
    print(f"FIBERPARAM {cart} {cam} {opFibers['mjd'][0]} 70.0000 {np.sum(opFibers['bundlefibers'][0])} "+
          f"{opFibers['nbundles'][0]} {fiberspacing} {bundlegaps} {bundlefibers} {deadfibermask}")
    return(xcen)
