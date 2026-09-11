#!/usr/bin/env python3
import numpy as np
from boss_drp.utils.merge_ranges import merge_ranges

def mjd_match(thismjd, mjd=None, mjdstart=None, mjdend=None, ranges=None):
    if mjd is not None:
        if int(thismjd) not in np.atleast_1d(np.asarray(mjd)).astype(int).tolist():
            return(False)
    if mjdstart is not None:
        if isinstance(mjdstart, (list,tuple)):
            valid = False
            for r in mjdstart:
                if int(thismjd) >= int(r[0]) and int(thismjd) <= int(r[1]):
                    valid = True
            if not valid:
                return(False)
        elif int(thismjd) < int(mjdstart):
                return(False)
    if mjdend is not None:
        if isinstance(mjdend, (list,tuple)):
            valid = False
            for r in mjdend:
                if int(thismjd) >= int(r[0]) and int(thismjd) <= int(r[1]):
                    valid = True
            if not valid:
                return(False)   
        elif int(thismjd) > int(mjdend):
            return(False)
    if ranges is not None:
        if isinstance(ranges[0], (list,tuple)):
            valid = False
            ranges = merge_ranges(ranges)
            for r in ranges:
                if r[0] is None:
                    r[0] = -9999
                if r[1] is None:
                    r[1] = 999999999
                if int(thismjd) >= int(r[0]) and int(thismjd) <= int(r[1]):
                    valid = True
            if not valid:
                return(False)
        else:
            if ranges[0] is None:
                ranges[0] =  -9999
            if ranges[1] is None:
                ranges[1] = 999999999
            if int(thismjd) < int(ranges[0]) or int(thismjd) > int(ranges[1]):
                return(False)
    return(True)
