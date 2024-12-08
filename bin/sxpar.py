#!/usr/bin/env python3
from boss_drp.utils.sxpar import sxpar
import argparse
import sys
import os


""" 
sxpar:

Simply parse a fits header.  Copied from perl "sxpar by D. Finkbeiner 2001 Dec 20".

Can read uncompressed or gz files.

Written by Gary Kushner (LBL).  Oct 2009.

"""

    
if __name__=='__main__':
    parser = argparse.ArgumentParser(
            prog=os.path.basename(sys.argv[0]),
            description='Simply parse a fits header')
    parser.add_argument('fitsfile',help='The fits file to read')
    parser.add_argument('keyword',help='Header keyword to parse')
    parser.add_argument('-v','--verbose', action='store_true', help='verbose')
    args = parser.parse_args()

    output = sxpar(args.fitsfile, args.keyword, args.verbose)
    for l in output:
        print(l)
    

