#!/usr/bin/env python3

from slurm import queue
import argparse
from os import getenv, makedirs, popen, chdir, getcwd
import os.path as ptt
from datetime import date
import astropy.time
import subprocess
import io
import sys
import pandas as pd
from load_module import load_module, load_env
from dailylogger import *
import logging
from pydl.pydlutils.yanny import yanny
from astropy.io import fits
from astropy.table import Table
import numpy as np
import time
import re


queue = queue()
mjd = str(int(float(astropy.time.Time(str(date.today())).jd) - 2400000.5))

#def load_env(key,log):
#    val = getenv(key)
#    if val is None:
#        log.info('ERROR: '+key+' is not set')
#        exit()
#    return(val)

def read_mod(mod,log):
    module = load_module()
    module('purge')
    module('load', mod)
    boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX',log)
    run2d = load_env('RUN2D',log)
    scratch_dir = load_env('SLURM_SCRATCH_DIR',log)
    return(boss_spectro_redux,run2d,scratch_dir)

def check_daily(mod, daily_dir, mjd, log):
    nextmjds = yanny(ptt.join(daily_dir, 'etc', 'nextmjd.par'))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where(mods == mod.lower())[0]
    if len(indx) == 0:
        return(False)
    else:
        return(nextmjds['NEXTMJD']['mjd'][indx].max() > mjd)

def check_fieldlist(boss_spectro_redux, run2d, spall_mjd):
    flist = Table(fits.getdata(ptt.join(boss_spectro_redux, run2d, 'fieldlist-'+run2d+'.fits'),1))
    r = re.compile('Done[\w]*', re.IGNORECASE)
    idx = [i for i, x in enumerate(flist['STATUS1D'].data ) if r.search(x)]
    flist = flist[idx]
    latest_redux = max(flist['MJD'])
    return(latest_redux > spall_mjd)

def slurm_fieldmerge(module='bhm/master', walltime = '40:00:00', fast = True, mem = 32000, daily=False):

    daily_dir = getenv('DAILY_DIR')
    if daily_dir is None: daily_dir = ptt.join(getenv('HOME'), "daily")

    makedirs(ptt.join(daily_dir, "logs", "fieldmerge", 'control'), exist_ok = True)
    filelog = logging.FileHandler(ptt.join(daily_dir, 'logs', 'fieldmerge', 'control', str(mjd)+'.log'))
    filelog.setLevel(logging.DEBUG)
    filelog.setFormatter(Formatter())

    logger = logging.getLogger()
    elog = emailLogger()
    emaillog = elog.log_handler
    emaillog.setLevel(logging.DEBUG)
    emaillog.setFormatter(Formatter())
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(Formatter())
    logger.addHandler(console)
    logger.addHandler(emaillog)
    logger.addHandler(filelog)
    logger.setLevel(logging.DEBUG)

    logger.info('===============================================')
    logger.debug(time.ctime())
    queue.verbose = True
    if fast is True: 
        alloc = 'sdss-np-fast'
    else:
        alloc = 'sdss-np'
    boss_spectro_redux, run2d, slurm_dir = read_mod(module, logger)

    if daily is True: 
        with fits.open(ptt.join(boss_spectro_redux, run2d, 'spAll-'+run2d+'.fits')) as ff:
            tf = Table(ff[1].data)
            latest_mjd = tf['MJD'].max()
        
        if not check_daily(module, daily_dir, latest_mjd, logger):
            logger.debug('Skipping run')
            elog.send('fieldmerge '+run2d +' MJD='+str(mjd), ptt.join(daily_dir, 'etc','emails'), logger)
            return()
        if not check_fieldlist(boss_spectro_redux, run2d, latest_mjd):
            logger.debug('SpAll-'+run2d+' up to date')
            elog.send('fieldmerge '+run2d +' MJD='+str(mjd), ptt.join(daily_dir, 'etc','emails'), logger)
            return()
    log = ptt.join(daily_dir, "logs", "fieldmerge", run2d, "fieldmerge_"+mjd)

    pdargs = pd.Series({"module":module,
                        "walltime":walltime,
                        "fast":fast,
                        "mem":mem,
                        "run2d":run2d})
    logger.info(pdargs.to_string())

    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout

    makedirs(ptt.join(daily_dir, "logs", "fieldmerge", run2d), exist_ok = True)
    title = 'fieldmerge_'+run2d
    queue.create(label = title, nodes = 1, ppn = 1, walltime = walltime, alloc=alloc, qos  = 'sdss',  partition = 'sdss-np',
                 mem_per_cpu = mem, shared = True)

    queue.append("module purge ; module load "+module+" ; source "+ptt.join(daily_dir, "cmd", 'run_fieldmerge'),
                 outfile = log+".o.log", errfile = log+".e.log")
    queue.commit(submit=True)

    output = new_stdout.getvalue()
    sys.stdout = old_stdout
    logger.info(output)
    
    subject = 'fieldmerge '+module +' MJD='+str(mjd)
    elog.send(subject, ptt.join(daily_dir, 'etc','emails'), logger)
    logger.removeHandler(emaillog)
    emaillog.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create daily field merge slurm job')
    parser.add_argument('--module', '-m', default = 'bhm/master', help = 'module file to use (ex bhm/master[default] or bhm/v6_0_9)')
    parser.add_argument('--walltime', '-w', default = '40:00:00', help = 'Job wall time (format hh:mm:ss) default = "40:00:00"')
    parser.add_argument('--fast', action='store_true' , help = 'use fast allocation')
    parser.add_argument('--mem', default = 32000, help = 'memory in bytes')
    parser.add_argument('--daily', action = 'store_true', help = 'only run if daily run has been run today')
    args = parser.parse_args()
    slurm_fieldmerge(module=args.module, walltime = args.walltime, fast = args.fast, mem = args.mem, daily = args.daily)
