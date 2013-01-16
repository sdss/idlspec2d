#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# $Id$
#
"""Finds files that should be present but are not.
"""
#
# Imports
#
#import numpy as np
import pyfits
import os
import os.path
#
#
#
def main():
    """Main program.
    """
    platesfile = os.path.join(os.getenv('SPECTRO_REDUX'),'plates-SDSS-dr9.fits')
    p = pyfits.open(platesfile)
    plates = p[1].data
    goodplates = plates['PLATEQUALITY'] == 'marginal'
    for run2d,plate,mjd in zip(plates['RUN2D'][goodplates],plates['PLATE'][goodplates],plates['MJD'][goodplates]):
        pmjd = "{0:04d}-{1:05d}".format(int(plate),int(mjd))
        files1 = [os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'{0}-{1}.fits'.format(f,pmjd)) for f in ('spPlate','spZall','spZbest','spZline','spFluxdistort')]
        for f in files1:
            if not os.path.isfile(f):
                print("{0} MISSING!".format(f))
        print("Examining {0}...".format(files1[0]))
        sp = pyfits.open(files1[0])
        hdr = sp[0].header
        expid = list()
        while hdr.has_key("EXPID{0:02d}".format(len(expid)+1)):
            foo = hdr["EXPID{0:02d}".format(len(expid)+1)]
            expid.append(foo)
        sp.close()
        for e in expid:
            ee = e.split('-')
            files2 = [os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'{0}-{1}-{2}.fits'.format(f,ee[0],ee[1])) for f in ('spFrame','spCFrame','spFluxcorr','spFluxcalib')]
            for f in files2:
                if not os.path.isfile(f) and not os.path.isfile(f+'.gz'):
                    print("{0} MISSING!".format(f))
            flat = os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'spFlat-{0}-{1}.fits.gz'.format(ee[0],ee[2]))
            if not os.path.isfile(flat):
                print("{0} MISSING!".format(flat))
            arc = os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'spArc-{0}-{1}.fits.gz'.format(ee[0],ee[3]))
            if not os.path.isfile(arc):
                print("{0} MISSING!".format(arc))
    p.close()
    return
#
#
#
if __name__ == '__main__':
    main()
