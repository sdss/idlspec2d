#!/usr/bin/env python

import os
import sys
from boss_drp.utils import sxpar
from boss_drp.utils import putils

"""
filecheck.py:

filecheck.py checks some attribute of a file and returns "true" or "false" based
on if the attribute is...  true or false.

This file it pretty much a grab bag of hacks.  To do these tests correctly would
either require adding dependencies just for these simple tests or duplicating lots
of code.  

The input files can be either normal or gz files.

Written by Gary Kushner (LBL).  Dec 2009.

"""

####
def science(fits, check_fix=False):
    """return True if the fits file is a science frame"""
    v = sxpar.sxparRetry(fits, "flavor", retries = 60, check_fix=check_fix)
    if len(v) == 0:
        return False
    return v[0].lower() == "science"
    


####
def excellent(fits, retries=60, return_qaulity=False, check_fix=False):
    
    """return True if the fits file is an excellent frame"""
    v = sxpar.sxparRetry(fits, "quality", retries = retries, check_fix=check_fix)
    if len(v) == 0:
        return True, 'excellent'
    if return_qaulity:
        return v[0].lower() == "excellent", v[0]
    return v[0].lower() == "excellent"

####
def arc(fits, check_fix=False):
    """return True if the fits file is a arc frame"""
    v = sxpar.sxparRetry(fits, "flavor", retries = 60, check_fix=check_fix)
    if len(v) == 0:
        return False
        
    if v[0].lower() == "arc":
        hart = sxpar.sxparRetry(fits, "HARTMANN", retries = 60, check_fix=check_fix)
        if len(hart) != 0:
            if hart[0].lower() in ['right','left']:
                return False
        
    return v[0].lower() == "arc"
####
def boss(fits):
    """return True if the plugmap file is for a boss frame"""
    
    key   = "instruments"
    keyl  = len(key)
    value = "boss"
    
    with putils.openRead(fits) as f:  # Ensure file is opened in the context
        for line in f:  # Read line by line (more memory efficient than readlines)
            if line[:keyl].lower() == key:
                return line.lower().count(value) > 0  # Returns immediately, but still closes file

    return False  # If no match is found
        
####
def test(fits, check_fix=False):
    """return True if the fits file is a test frame"""

    v = sxpar.sxparRetry(fits, "quality", retries = 60, check_fix=check_fix)
    if len(v) == 0:
        return False
    return v[0].lower() == "test"


def filecheck(cmd, file):
    if cmd == "science":
        rv = science(file)
    elif cmd == "arc":
        rv = arc(file)
    elif cmd == "test":
        rv = test(file)
    elif cmd == "excellent":
        rv = excellent(file)
    elif cmd == "boss":
        rv = boss(file)
    else:
        rv = False
    if rv == True:
        print("true")
    else:
        print("false")
