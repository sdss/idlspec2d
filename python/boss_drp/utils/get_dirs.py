#!/usr/bin/env python3
from boss_drp.field import field_to_string
from boss_drp.utils.merge_ranges import merge_ranges

import os.path as ptt
import numpy as np
import glob


def get_dirs(basedir, subdir='', pattern='*', match=None,
            start=None, end=None, ranges=None, numeric=True, field=False):
    """
        Generates of list of directores matching a patten with in basedir/subdir, and filters out
        folders outside of valid range
    """
    dlist = []
    dirs =  glob.glob(ptt.join(basedir, subdir, pattern))

    if field is True:
        if pattern == '*':
            pathern = field_to_string('*')
    if match is not None:
        if numeric:
            match = np.atleast_1d(np.asarray(match)).astype(int).tolist()
        else:
            match = np.atleast_1d(np.asarray(match)).astype(str).tolist()
    
    for d in dirs:
        d = ptt.basename(d)
        if numeric is True:
            if not ptt.basename(d).isnumeric():
                continue
        if match is not None:
            if numeric:
                if int(d) not in match:
                    continue
            else:
                if d not in match:
                    continue
        if ranges is not None:
            if isinstance(ranges[0], (list,tuple)):
                ranges = merge_ranges(ranges)
                valid = False
                for r in ranges:
                    if r[0] is None:
                        r[0] = -10000000
                    if r[1] is None:
                        r[1] = 9999999999
                    if int(d) >= int(r[0]) and int(d) <= int(r[1]):
                        valid = True
                if not valid:
                    continue
            else:
                if int(d) < int(ranges[0]) or int(d) > int(ranges[1]):
                    continue
        if start is not None:
            if int(d) < int(start):
                continue
        if end is not None:
            if int(d) > int(end):
                continue
        dlist.append(d)
    dlist = np.sort(np.asarray(dlist)).tolist()
    return(dlist)

