#!/usr/bin/env python3

from boss_drp.summary import summary_names as fnames
from boss_drp.utils.splog import splog
import os
import numpy as np
from collections import OrderedDict
from glob import glob
from datetime import datetime

def cleanup_bkups(boss_spectro_redux, run2d, backups = 3, epoch=False,
                  custom = None):
    allsky = False if custom is not None else True
    fnames.set(boss_spectro_redux, run2d, dev=False,
                           epoch=epoch, allsky=allsky, custom=custom)
    bk_files = OrderedDict()
    fnames.bk.set(flag = '*')
    
    # Gather backup files and extract timestamp
    for bf in glob(fnames.bk.spAllfile):
        try:
            timestamp_str = bf.split('.bkup-')[-1]
            timestamp = datetime.strptime(timestamp_str, "%Y-%m-%dT%H:%M:%S")
            bk_files[bf] = timestamp
        except Exception as e:
            splog.warning(f"Skipping file {bf}: could not parse timestamp ({e})")
   
    # Handle no backups
    if len(bk_files) == 0:
        splog.debug('No Backups exist')
        return

    # Sort by datetime, descending (newest first)
    sorted_keys = list(bk_files.keys())
    sorted_timestamps = list(bk_files.values())
    idx = np.argsort(sorted_timestamps)[::-1]  # Reverse for newest first

    # Cleanup logic
    for i, id in enumerate(idx):
        key = sorted_keys[id]
        if i > backups - 1:
            splog.debug(f'Removing Backup: {key} ({i+1})')
            os.remove(key)
            for variant in ['spAll-lite-', 'spAllLine-']:
                alt_file = key.replace('spAll-', variant)
                splog.debug(' '*16 + alt_file)
                os.remove(alt_file)
        else:
            splog.debug(f'Keeping Backup: {key} ({i+1})')

    return
