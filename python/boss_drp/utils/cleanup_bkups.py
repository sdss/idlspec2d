#!/usr/bin/env python3

from boss_drp.summary import summary_names as fnames
from boss_drp.utils.splog import splog
import os
import numpy as np
from collections import OrderedDict
from glob import glob

def cleanup_bkups(boss_spectro_redux, run2d, backups = 3, epoch=False,
                  custom = None):
    allsky = False if custom is not None else True
    fnames.set(boss_spectro_redux, run2d, dev=False,
                           epoch=epoch, allsky=allsky, custom=custom)
    bk_files = OrderedDict()
    fnames.bk.set(flag = '*')
    for bf in glob(fnames.bk.spAllfile):
        bk_files[bf] = bf.split('-')[-1]
    if len(bk_files) == 0:
        splog.debug('No Backups exist')
    idx = np.argsort(np.asarray(list(bk_files.values())))
    for i, id in enumerate(np.flip(idx)):
        key = list(bk_files.keys())[id]
        if i > backups -1:
            splog.debug(f'Removing Backup: {key} ({i+1})')
            os.remove(key)
            splog.debug(' '*16+key.replace('spAll-','spAll-lite-'))
            os.remove(key.replace('spAll-','spAll-lite-'))
            splog.debug(' '*16+key.replace('spAll-','spAllLine-'))
            os.remove(key.replace('spAll-','spAllLine-'))
        else:
            splog.debug(f'Keeping Backup: {key} ({i+1})')
    return
