#!/usr/bin/env python3

import sys
import os
import time
from boss_drp.utils import putils
from boss_drp.utils import Sphdrfix
"""
sxpar:

Simply parse a fits header.  Copied from perl "sxpar by D. Finkbeiner 2001 Dec 20".

Can read uncompressed or gz files.

Written by Gary Kushner (LBL).  Oct 2009.

"""
    
####
def sxparRetry(fitsfile, keyword = None, verbose = False, retries = 0, check_fix=False):
    """call sxpar with retries.  To see if the fits wasn't finished writing."""
    while True:
        try:
            return sxpar(fitsfile, keyword, verbose, check_fix=check_fix)
        except:
            if retries < 0:
                raise
            retries = retries - 1
            time.sleep(1)
            
####
def sxpar(fitsfile, keyword = None, verbose = False, mjd=None, check_fix=False):
    """Parse fits header and return output list"""
    
    isFits = False
    output = []

    if keyword != None:
        keyword = keyword.upper()
    
    if (mjd is None) and (check_fix):
        mjd = sxpar(fitsfile, keyword='MJD')
        obs = sxpar(fitsfile, keyword='CARTID')
        obs = 'LCO' if str(obs).strip().lower() == 'fps-s' else 'APO'
        hdrfix = Sphdrfix(mjd, obs = obs)
        
    with putils.openRead(fitsfile, mode='rb') as f:
        i = 0
        while True:
            i += 1
            line = f.read(80).decode("utf-8")
            
            key = line.split("=")[0].strip()
            
            if key == "SIMPLE":
                isFits = True
                continue
            if key == "END":
                break

            if i > 40 and not isFits:
                raise TypeError(fitsfile + " Doesn't look like a fits file -- did not find 'SIMPLE'")

            if keyword == None or key == keyword:
                values  = line.partition("=")[2].partition("/")
                value   = values[0].strip().strip("'").strip()
                if check_fix:
                    try:
                        _fix = hdrfix.fix(fitsfile, {key:value})
                        value = _fix[key]
                    except:
                        pass
                if verbose:
                    output.append(line)
                elif keyword == None:
                    output.append(key + " = " + value)
                else:
                    output.append(value)
    
    return output
