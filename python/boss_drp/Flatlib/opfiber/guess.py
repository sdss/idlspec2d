from boss_drp import idlspec2d_dir

from astropy.io import fits
from scipy.signal import find_peaks
from pydl.pydlutils.yanny import read_table_yanny

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os.path as ptt
import warnings

def build_trace_guess(procFits, bundlefibers = None, mjd = None, plot=False, min_peak_sep = 6, min_peak_height=5000, precision=3):
    """ Takes a processed image frame produced using the /sawraw flag in sdssproc (or indirecly via spreduce2d)
        as an input. It then uses either the bundlefiber list of number of fibers per bundle (supplied as input or via opFiberFPS)
        combined with the scipy peak finding algarithm to create a first guess of the peak fiber and bundle gaps. It uses the median
        flux of the 11 central pixel (along the dispersion axis) to build the flux array
        
        procFits : str
            Flat fits file name of form sdProc-XX-XXXXXXXX.fits
        bundlefibers : list(int) (optional
            list of number of fibers per bundle (defaults to reading from latest opFiberFPS entry)
        mjd : int (optional)
            MJD of new updated opFiber Fiberparameter entry
        plot : bool
            Whether to plot the flux slice and detected peaks
        min_peak_sep : float (optional: default=6)
            the minimum seperation between detected peak (ie fibers are no closer then this)
        min_peak_height : float (optional : default=5000)
            the minimum flux level to be detected as a peak
            
    """
    
    matplotlib.use('Agg')  # Use non-interactive backend

    flux = fits.getdata(procFits,0)
    ny = flux.shape[0]
    ystart = ny//2
    nmed = 11
    ylo = max(0, int(ystart - (nmed - 1) / 2))
    yhi = min(ny - 1, int(ystart + (nmed - 1) / 2))

    flux = np.median(flux[ylo:yhi+1,:],axis=0)
    xpix = np.arange(len(flux))  # Unshifted x-axis for fluxvec
    peaks,_ = find_peaks(flux, distance=min_peak_sep, height=min_peak_height)
    xcen = xpix[peaks]

    if plot:
        nplotrow = 8
        fig, axs = plt.subplots(nplotrow, 1, figsize=(8,11))#, sharex=True)
        for j in range(nplotrow):
        # Define consistent xrange (before shifting)
            xrange = len(flux) * np.array([j, j + 1]) / float(nplotrow)
            xrange[1] += 5  # Slightly extend the upper bound
            axs[j].plot(xpix, flux, color='k', label = ptt.basename(procFits).replace('.fits','').replace('sdProc-',''))
            axs[j].plot(xcen, flux[peaks], ls = '', marker ='x', color='r', label='peak')
            axs[j].set_xlim(xrange)
        axs[0].legend()
        plt.tight_layout()
        plt.savefig(ptt.basename(procFits.replace('.fits', '.jpeg').replace('sdProc-', 'spPeaks-')))
        plt.show()

    cam  = ptt.basename(procFits).split('-')[1]
    cart = 'FPS-S' if '2' in cam else 'FPS-N'

    
    opFibers = read_table_yanny(ptt.join(idlspec2d_dir, 'opfiles', 'opFibersFPS.par'), 'FIBERPARAM')
    opFibers = opFibers[(opFibers['cartid'] == cart) & (opFibers['camname'] ==cam)]
    opFibers= opFibers[opFibers['mjd'] == opFibers['mjd'].max()]
    opFibers
    
    bundlegaps = []
    fiberspacing = []
    sfiber = 0
    if bundlefibers is None:
        bundlefibers = opFibers['bundlefibers'][0]
    for i, nfib in enumerate(bundlefibers):
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
    
    bundlefibers_s = '{ '+" ".join(f"{x:d}" for x in bundlefibers)+' }'
    deadfibermask = '{ '+" ".join(f"{x:d}" for x in opFibers['deadfibermask'][0])+' }'

    print('--------------------- New version ----------------------')

    if mjd is None:
        mjd = opFibers['mjd'][0]
    nbund = np.count_nonzero(bundlefibers)
    nfibers = np.sum(bundlefibers)
    print(f"FIBERPARAM {cart} {cam} {mjd} 70.0000 {nfibers} "+
          f"{nbund} {fiberspacing} {bundlegaps} {bundlefibers_s} {deadfibermask}")
    
    print('------------------ Previous version --------------------')
    fiberspacing = '{ '+" ".join(f"{x:{format_str}}" for x in opFibers['fiberspace'][0])+' }'
    bundlegaps = '{ '+" ".join(f"{x:{format_str}}" for x in opFibers['bundlegap'][0])+' }'
    bundlefibers = '{ '+" ".join(f"{x:d}" for x in opFibers['bundlefibers'][0])+' }'
    print(f"FIBERPARAM {cart} {cam} {opFibers['mjd'][0]} 70.0000 {np.sum(opFibers['bundlefibers'][0])} "+
          f"{opFibers['nbundles'][0]} {fiberspacing} {bundlegaps} {bundlefibers} {deadfibermask}")
    return
