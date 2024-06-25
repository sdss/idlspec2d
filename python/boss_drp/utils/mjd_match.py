#!/usr/bin/env python3
import numpy as np

def mjd_match(thismjd, mjd=None, mjdstart=None, mjdend=None):
    if mjd is not None:
        if int(thismjd) not in np.atleast_1d(np.asarray(mjd)).astype(int).tolist():
            return(False)
    if mjdstart is not None:
        if int(thismjd) < int(mjdstart):
            return(False)
    if mjdend is not None:
        if int(thismjd) > int(mjdend):
            return(False)
    return(True)
