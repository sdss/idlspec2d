#!/usr/bin/env python

"""
Try solving with a prior that fluxcorr = 1
"""

import sys
import os
import os.path
import numpy as N
import pylab as P
import fitsio
from numpy.polynomial import chebyshev
from scipy.sparse.construct import spdiags
import yanny

#-------------------------------------------------------------------------

def read_plan(planfile):
    plan = yanny.read_yanny(planfile)
    plandir = os.path.dirname(planfile)
    framefiles = dict(b1=list(), b2=list(), r1=list(), r2=list())

    for i in range(len(plan['SPEXP']['name'])):
        rawfiles = plan['SPEXP']['name'][i]
        for filename in rawfiles:
            if 'UNKNOWN' in filename: continue
            pre, camera, exp = os.path.splitext(filename)[0].split('-')
            if not os.path.exists(filename) and not filename.endswith('.gz'):
                filename = filename + '.gz'
            if not os.path.exists(filename):
                print 'fluxcorr_prior.py: File does not exist, skipping. %s'%filename
                continue
            framefiles[camera].append(filename)

    return framefiles

def read_data(indir, framefiles, xythrucorr=False):
    flux = list()
    ivar = list()
    goodframes = set()
    
    for framefile in framefiles:
        infile = indir + '/' + framefile
            
        calibfile = infile.replace('spFrame', 'spFluxcalib')
        xythrufile = infile.replace('spFrame', 'spXYthrucorr')

        if not os.path.exists(infile):
            print "WARNING: Skipping missing", infile
            continue

        if not os.path.exists(calibfile):
            print "WARNING: Skipping missing", calibfile
            continue

        if xythrucorr and not os.path.exists(xythrufile):
            print "WARNING: Skipping missing", xythrufile
            continue

        goodframes.add(framefile)

        calib = fitsio.read(calibfile, 0)
        goodcalib = (calib != 0)

        print "Reading", os.path.basename(infile)
        fx = fitsio.FITS(infile)
        xflux = fx[0].read()
        xflux[goodcalib] /= calib[goodcalib]

        xivar = fx[1].read()
        xmask = fx[2].read()
        xivar[xmask != 0] = 0.0
        xivar[goodcalib] *= calib[goodcalib]**2

        if xythrucorr:
            print "Reading", os.path.basename(xythrufile)
            thrucorr = fitsio.read(xythrufile, 0)
            xflux[goodcalib] *= thrucorr[goodcalib]
            xivar[goodcalib] /= (thrucorr[goodcalib])**2


        #- Add to flux and ivar lists
        flux.append(xflux)
        ivar.append(xivar)

        fx.close()

    flux = N.array(flux)
    ivar = N.array(ivar)
    flux[ivar==0] = 0.0  #- cosmetics

    return flux, ivar, goodframes
    
    
def calc_fluxcorr(flux, ivar, prior=0.5):

    print "Calculating flux corrections"
    
    #- The array to fill
    fluxcorr = N.ones(flux.shape)  
    
    #- Determine range to actually fill
    ii = N.where( N.sum(N.sum(ivar, axis=0), axis=0) > 0 )[0]
    imin, imax = ii[0], ii[-1]+1
    flux = flux[:, :, imin:imax]
    ivar = ivar[:, :, imin:imax]
    
    nexp, nspec, npix = flux.shape
    
    #- Make coadd; handle cases where sum(weights) = 0
    weighted_flux = N.sum(flux*ivar, axis=0)
    sum_weights = N.sum(ivar, axis=0)
    coadd = N.zeros(weighted_flux.shape)
    ii = (sum_weights > 0)
    coadd[ii] = weighted_flux[ii] / sum_weights[ii]

    #- Create Chebyshev matrix
    #- FL = diag(coadd).dot(ChebyPolys)
    #- flux[i] = FL.dot(c[i])
    npoly = 4
    xx = N.linspace(-1, 1, npix)

    L = N.zeros( (npix, npoly) )
    for i in range(npoly):
        c = N.zeros(npoly)
        c[i] = 1.0
        L[:,i] = chebyshev.chebval(xx, c)
    
    for ispec in range(nspec):
        c = list()
        FL = spdiags(coadd[ispec], 0, npix, npix).dot(L)
        for iexp in range(nexp):
            #- Create diagonal weights matrix, throwing out worst 5% for robustness
            weights = ivar[iexp, ispec].copy()
            if N.sum(weights) == 0:
                corr = N.zeros(npoly)
                corr[0] = 1.0
                c.append(corr)
                continue

            wcut = N.percentile(weights[weights>0], 5)        
            weights[weights<wcut] = 0.0
            Wi = spdiags(weights, 0, npix, npix)
            
            #- Weighted flux
            wf = Wi.dot(flux[iexp, ispec])
    
            #- Weighted coadd*chebyshev polynomials
            WFL = Wi.dot(FL)
    
            #- Original constant prior that fluxcorr=1
            # wf = N.concatenate( (wf, prior*N.ones(npix)) )
            # WFL = N.vstack( (WFL, prior*L) )
            
            #- Add prior that fluxcorr = 1, scaled by data weights
            wf = N.concatenate( (wf, prior*weights) )            
            diagtmp = spdiags(prior*weights, 0, npix, npix)
            WFL = N.vstack( (WFL, diagtmp.dot(L)) )
     
            #- Solve for chebyshev coefficients
            c.append( N.linalg.lstsq(WFL, wf)[0] )
            
            fluxcorr[iexp, ispec, imin:imax] = 1.0 / chebyshev.chebval(xx, c[iexp])

    return fluxcorr

def plotstuff(flux, fluxcorr, ispec):
    nexp, nspec, npix = flux.shape
    
    xx = N.linspace(-1, 1, npix)
    ii = N.where(N.sum(N.sum(flux, axis=0), axis=0) != 0)[0]
    xmin, xmax = xx[ii[0]], xx[ii[-1]]
    
    P.clf()
    P.subplot(311)
    for iexp in range(nexp):
        P.plot(xx, 1/fluxcorr[iexp, ispec])
    P.ylim(0,2)
    P.xlim(xmin, xmax)

    P.subplot(312)
    window = N.hamming(51)
    window /= N.sum(window)
    for iexp in range(nexp):
        P.plot(xx, N.convolve(flux[iexp, ispec], window, mode='same'))
    P.xlim(xmin, xmax)
        
    ymax = N.percentile(flux[:, ispec], 95)
    P.ylim(-2,ymax)
    P.xlim(xmin, xmax)

    P.subplot(313)
    for iexp in range(nexp):
        P.plot(xx, N.convolve(flux[iexp, ispec]*fluxcorr[iexp, ispec], window, mode='same'))
    P.ylim(-2,ymax)
    P.xlim(xmin, xmax)
    print ispec
    
#-------------------------------------------------------------------------
planfile = sys.argv[1]
if len(sys.argv) > 2:
	xythrucorr = True
else:
	xythrucorr = False
 
plandir = os.path.dirname(os.path.abspath(planfile))
framefiles = read_plan(planfile)
for camera in ('b1', 'b2', 'r1', 'r2'):
    flux, ivar, goodframes = read_data(plandir, framefiles[camera], xythrucorr=xythrucorr)
    if len(goodframes)>0:
        fluxcorr = calc_fluxcorr(flux, ivar, prior=1.0)
        addterm = N.zeros(fluxcorr[0].shape)
   
    i = 0
    for framefile in framefiles[camera]:
        corrfile = plandir+'/'+framefile.replace('spFrame', 'spFluxcorr')
        if corrfile.endswith('.gz'):
            corrfile = corrfile[:-3]
           
        print "Writing", os.path.basename(corrfile)
        if framefile in goodframes:
            corr=fluxcorr[i]
            add = addterm
            i+=1
        else:
            fluxshape = fitsio.read(framefile+'.gz', 0).shape
            corr = N.ones(fluxshape)
            add = N.zeros(fluxshape)

        fitsio.write(corrfile, corr, clobber=True)
        fitsio.write(corrfile, add)
        os.system('gzip -f '+corrfile)

# P.ion()
# i = 0
# plotstuff(flux, fluxcorr, i)

