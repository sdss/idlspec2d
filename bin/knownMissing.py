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
from __future__ import print_function
#import numpy as np
import pyfits
import argparse
import os
import os.path
import sys
#
# Top-level definitions
#
__author__ = 'Benjamin Alan Weaver <benjamin.weaver@nyu.edu>'
__version__ = '$Revision$'.split(': ')[1].split()[0]
__docformat__ = 'restructuredtext en'
#
# Functions
#
def version(headurl="$HeadURL$"):
    """Returns a version string based on the value of the svn HeadURL keyword.

    Parameters
    ----------
    headurl : str, optional
        A HeadURL string set by svn.

    Returns
    -------
    myversion : str
        A version string.

    Examples
    --------
    >>> transfer.common.version()
    trunk
    """
    if headurl.find('$HeadURL:') < 0:
        raise ValueError("The input does not appear to be a HeadURL string!")
    if headurl.find('trunk') > 0:
        myversion = 'trunk'
    elif headurl.find('tags') > 0:
        myversion = headurl[headurl.find('tags')+5:].split('/')[0]
    elif headurl.find('branches') > 0:
        myversion = headurl[headurl.find('branches')+9:].split('/')[0]
    else:
        myversion = 'unknown'
    return myversion
#
#
#
def main():
    """Main program.
    """
    #
    # Get options
    #
    parser = argparse.ArgumentParser(description=__doc__,prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-V','--version',action='version',
        version=("%(prog)s {0}".format(version())),
        help='Print version information and exit.')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
        help='Print extra information.')
    #parser.add_argument('-t', '--test', action='store_true', dest='test',
    #    help='Do not actually run commands (implies --verbose).')
    parser.add_argument('-m', '--marginal', action='store_true', dest='marginal',
        help='Analyze marginal plates instead of good plates.')
    parser.add_argument('-r', '--release', action='store', dest='release',
        default=9, type=int, help='Assume this data release.',
        metavar='RELEASE')
    options = parser.parse_args()
    #
    # Read the plates file.
    #
    platesfile = os.path.join(os.getenv('SPECTRO_REDUX'),'plates-SDSS-dr{0:d}.fits'.format(options.release))
    p = pyfits.open(platesfile)
    plates = p[1].data
    if options.marginal:
        goodplates = plates['PLATEQUALITY'] == 'marginal'
    else:
        goodplates = plates['PLATEQUALITY'] == 'good'
    for run2d,plate,mjd in zip(plates['RUN2D'][goodplates],plates['PLATE'][goodplates],plates['MJD'][goodplates]):
        pmjd = "{0:04d}-{1:05d}".format(int(plate),int(mjd))
        files1 = [os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'{0}-{1}.fits'.format(f,pmjd)) for f in ('spFluxdistort','spPlate','spZall','spZbest','spZline')]
        for f in files1:
            if not os.path.isfile(f):
                print("{0} MISSING!".format(f))
        if options.verbose:
            print("Examining {0}...".format(files1[1]))
        sp = pyfits.open(files1[1])
        hdr = sp[0].header
        expid = list()
        while hdr.has_key("EXPID{0:02d}".format(len(expid)+1)):
            foo = hdr["EXPID{0:02d}".format(len(expid)+1)]
            expid.append(foo)
        sp.close()
        for e in expid:
            ee = e.split('-')
            arc = os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'spArc-{0}-{1}.fits.gz'.format(ee[0],ee[3]))
            if not os.path.isfile(arc):
                print("{0} MISSING!".format(arc))
            cframe = os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'spCFrame-{0}-{1}.fits'.format(ee[0],ee[1]))
            if not os.path.isfile(cframe):
                print("{0} MISSING!".format(cframe))
            flat = os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'spFlat-{0}-{1}.fits.gz'.format(ee[0],ee[2]))
            if not os.path.isfile(flat):
                print("{0} MISSING!".format(flat))
            files2 = [os.path.join(os.getenv('SPECTRO_REDUX'),run2d,'{0:04d}'.format(int(plate)),'{0}-{1}-{2}.fits.gz'.format(f,ee[0],ee[1])) for f in ('spFluxcalib','spFluxcorr','spFrame')]
            for f in files2:
                if not os.path.isfile(f):
                    print("{0} MISSING!".format(f))
    p.close()
    return
#
#
#
if __name__ == '__main__':
    main()
