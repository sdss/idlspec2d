#!/usr/bin/env python3
from boss_drp.utils import load_env
from boss_drp.utils import jdate
from boss_drp import daily_dir
from boss_drp.Config import config, Config
import sys
from boss_drp.run.queue import Queue

from os import getenv, makedirs, popen, chdir, getcwd
import os.path as ptt
from datetime import date, datetime
import astropy.time
import subprocess
import io
import numpy as np
from glob import glob

def set_flags(self, no_reject = False, clobber_fibermap = False,
                sdssv_sn2 = False, no_arc2trace = False, forcea2t = False,
                bright = False, sn2_15 = False):
    self.flags = ['--utah', '--nodb']
    if no_reject:
        self.flags.append('--no_reject')
    if clobber_fibermap:
        self.flags.append('-f')
    if sdssv_sn2:
        self.flags.append('--sdssv_sn2')
    if no_arc2trace:
        self.flags.append('-a')
    if forcea2t:
        self.flags.append('-o')
    if sn2_15:
        self.flags.append('--sn2_15')
    if bright:
        self.flags.append('--bright')

Config.flags = []
Config.set_flags = set_flags

def slurm_SOS(mjd=None, mjdstart=None, mjdend=None, obs=['apo','lco'], **flags):    
    qu = build(daily_dir=daily_dir, mjd=mjd, mjdstart= mjdstart,
               mjdend=mjdend, obs = obs, **flags)

def build(mjd=None, mjdstart= None, mjdend=None,
            obs = ['apo','lco'], daily_dir=daily_dir, **flags):
    i = 0
    title = 'SOS'
    log = ptt.join(daily_dir, "logs", "sos")
    makedirs(log, exist_ok = True)

    if mjdstart is not None:
        if mjdend is None:
            mjdend = jdate.astype(int)
        mjd = np.arange(mjdstart, mjdend+1)

    if mjd is None:
        mjd = jdate.astype(int)

    maxppn = len(mjd)*len(obs) #*len(['b','r'])
    if maxppn < config.queue.get('nodes')*config.queue.get('ntasks'):
        config.queue.set('ppn',maxppn)
        config.queue.set('cpus', 2 * maxppn)

    if not config.queue.get('no_submit'):
        print(config)
    config.set_flags(**flags)

    for mj in mjd:
        for ob in obs:
            for cam in ['j']:    #'b','r']:
                if i == 0:
                    queue1 = Queue(config.queue, key=None, verbose=True)
                    queue1.create(**config.queue.to_dict(label=title))

                thislog = ptt.join(log, f'SOS-{ob}-{cam}-{mj}')
                thiscmd = (f"export OBSERVATORY={ob.upper()} ; " +
                           f"SOS -{cam} -m {mj} "+" ".join(config.flags))

                queue1.append(thiscmd, outfile = thislog+".o.log", errfile = thislog+".e.log")
                i = i+1
    queue1.commit(hard=True,submit=config.queue.get('no_submit'))
    return(queue1)
