#!/usr/bin/env python3
from boss_drp.utils.sxpar import sxparRetry

import argparse
import sys
import os
"""
sxpar:

Simply parse a fits header.  Copied from perl "sxpar by D. Finkbeiner 2001 Dec 20".


"""

        
if __name__=='__main__':
    parser = argparse.ArgumentParser(
            prog=os.path.basename(sys.argv[0]),
            description='Simply parse a fits header, retrying if failed')
    parser.add_argument('fitsfile',help='The fits file to read')
    parser.add_argument('keyword',help='Header keyword to parse')
    parser.add_argument('-v','--verbose', action='store_true', help='verbose')
    args = parser.parse_args()

    output = sxparRetry(args.fitsfile, args.keyword, args.verbose, retry = 60)
    for l in output:
        print(l)
    
