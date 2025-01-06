#!/usr/bin/env python3
from boss_drp.post.fieldmerge import summary_names as fnames, fieldlist_name
from boss_drp.utils import jdate, send_email
from boss_drp.run import monitor_job
from boss_drp import daily_dir
from boss_drp.utils.splog import splog

try:
    from slurm import queue
    noslurm = False

except:
    import warnings
    class SlurmWarning(Warning):
        def __init__(self, message):
            self.message = message
    def __str__(self):
            return repr(self.message)
    warnings.warn('No slurm package installed: printing command to STDOUT for manual run',SlurmWarning)
    noslurm = True
    queue = None
    
from pydl.pydlutils.yanny import yanny
from astropy.io import fits
import astropy.time
from astropy.table import Table

import os
import os.path as ptt
from datetime import date
import io
import sys
import numpy as np
import time
import re
from glob import glob
from collections import OrderedDict


def check_daily(mod, daily_dir, mjd):
    nextmjds = yanny(ptt.join(daily_dir, 'etc', 'nextmjd.par'))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where(mods == mod.lower())[0]
    if len(indx) == 0:
        return(False)
    else:
        return(nextmjds['NEXTMJD']['mjd'][indx].max() > mjd)

def check_fieldlist(boss_spectro_redux, run2d, spall_mjd):
    fieldlist_name.build(boss_spectro_redux, run2d,epoch=False, custom_name=None)
    flist_file = fieldlist_name.name
    flist = Table(fits.getdata(flist_file,1))
    r = re.compile('Done[\w]*', re.IGNORECASE)
    idx = [i for i, x in enumerate(flist['STATUS1D'].data ) if r.search(x)]
    flist = flist[idx]
    latest_redux = max(flist['MJD'])
    return(latest_redux > spall_mjd)


class Setup:
    def __init__(self):
        self.module = None
        self.boss_spectro_redux = None
        self.run2d = None
        self.run1d = None
        self.alloc = None
        self.nodes = 1
        self.ppn = None
        #self.cpus = 1
        self.mem_per_cpu = None
        self.walltime = None
        self.shared = False
        self.partition = None
        self.epoch = False
        self.custom = None
        self.daily = False
        self.merge_only = False
        self.ndays = None
        self.backup = None
        self.limit = None
        self.n_iter = None
        self.no_fieldlist = False
        self.skip_specprimary = False
        self.update_specprimary = False
        self.verbose = False
        self.utah_daily = False

    def __repr__(self):
        return self.__str__()
                                                                                                        
    def __str__(self):
        return (f"boss_spectro_redux: {self.boss_spectro_redux} \n"    +
                f"run2d: {self.run2d} \n"    +
                f"run1d: {self.run1d} \n"    +
                f"epoch: {self.epoch} \n"    +
                f"custom: {self.custom} \n"    +
                f"alloc: {self.alloc} \n"    +
                f"partition: {self.partition} \n"    +
                f"nodes: {self.nodes} \n"    +
                f"ppn: {self.ppn} \n"    +
                #f"cpus: {self.cpus} \n"    +
                f"mem_per_cpu: {self.mem_per_cpu} \n"    +
                f"walltime: {self.walltime} \n"+
                f"shared: {self.shared} \n" +
                f"merge_only: {self.merge_only} \n" +
                f"ndays: {self.ndays} \n" +
                f"backup: {self.backup} \n" +
                f"limit: {self.limit} \n" +
                f"n_iter: {self.n_iter} \n" +
                f"no_fieldlist: {self.no_fieldlist} \n" +
                f"skip_specprimary: {self.skip_specprimary} \n" +
                f"update_specprimary: {self.update_specprimary} \n" +
                f"verbose: {self.verbose} \n" +
                f"daily: {self.daily} \n"+
                f"utah_daily: {self.utah_daily} \n")


def slurm_Summary(topdir, run2d, run1d = None, module = None, alloc=None, partition=None,
                  walltime = '40:00:00', fast = False, mem = None, daily=False,
                  epoch=False, custom=None, full=False, monitor=False, no_submit = False,
                  merge_only=True, backup=None, limit=None, n_iter=None,
                  email_start = False, no_fieldlist=False, skip_specprimary=False,
                  update_specprimary = False, verbose = False, ndays = None,
                  utah_daily = False):

    setup = Setup()
    setup.module = module
    if topdir is None:
        setup.boss_spectro_redux = os.getenv('BOSS_SPECTRO_REDUX')
    else:
        setup.boss_spectro_redux = topdir
    if run2d is None:
        setup.run2d = os.getenv('RUN2D')
    else:
        setup.run2d = run2d
    if run1d is None:
        setup.run1d = setup.run2d
    else:
        setup.run1d = run1d
    setup.daily = daily
    setup.epoch = epoch
    setup.custom = custom
    setup.merge_only = merge_only
    setup.backup = backup
    setup.limit = limit
    setup.n_iter = n_iter
    setup.no_fieldlist = no_fieldlist
    setup.skip_specprimary = skip_specprimary
    setup.update_specprimary = update_specprimary
    setup.verbose = verbose
    setup.ndays = ndays
    setup.utah_daily = utah_daily

    setup.alloc = alloc
    if setup.alloc is None:
        setup.alloc = os.getenv('SLURM_ALLOC')
    setup.partition = partition
    if setup.partition is None:
        setup.partition = setup.alloc
    setup.mem_per_cpu = os.getenv('SLURM_MEM_PER_CPU')
    if 'sdss-np' in setup.alloc:
        if fast is True:
            setup.alloc = 'sdss-np-fast'
        setup.shared = True
    else:
        full = True
        setup.shared=False
    
    if full:
        setup.shared = False
        setup.ppn = os.getenv('SLURM_PPN')
    else:
        setup.ppn = 10
        if mem is not None:
            setup.mem_per_cpu  = mem/setup.ppn
    setup.walltime = walltime
    
    queue1, title, attachements = build(setup, no_submit=no_submit,
                                        email_start = email_start)
                                
                                
    if no_submit:
        monitor = False
    if monitor:
        monitor_job(queue1, pause = 300, jobname = title)
        subject = _build_subject(setup, jdate.astype(str))
        #if not no_submit:
        #    cleanup_bkups(setup)
        splog.close_file()
        flags = []
        if setup.no_fieldlist:
            flags.append('Complete: fieldlist')
        if attachements is not None:
            for att in attachements:
                with open(att) as f:
                    lines = f.readlines()
                    lines.reverse()
                    for line in lines:
                        if 'Successful completion of fieldmerge' in line:
                            flags.append('Complete: fieldmerge')
                        elif 'Successful completion of fieldlist' in line:
                            flags.append('Complete: fieldlist')
                        elif 'exception:' in line:
                            flags.append('Crashed')
                        elif 'EXITING' in line:
                            flags.append('Killed')
                        elif 'Killed' in line:
                            flags.append('Killed')
            if 'Killed' in flags:
                subject = 'Killed: '+subject
            elif 'Crashed' in flags:
                subject = 'Crashed: '+subject
            elif ('Complete: fieldmerge' in flags and 'Complete: fieldlist' in flags):
                subject = 'Complete: '+subject
            else:
                subject = '???: '+subject
        send_email(subject, ptt.join(daily_dir, 'etc','emails'),
                      attachements, from_domain="chpc.utah.edu")

def _build_subject(setup, mjd):
    mstr = setup.module if setup.module is not None else setup.run2d
    if setup.epoch:
        subject = 'BOSS Summary '+mstr +'epoch MJD='+str(mjd)
    elif setup.custom is not None:
        subject = 'BOSS Summary '+mstr +f' {setup.custom} MJD='+str(mjd)
    else:
        subject = 'BOSS Summary '+mstr +' MJD='+str(mjd)
    return(subject)

def _build_log_dir(setup, control = False):
    daily_dir = os.getenv('DAILY_DIR')
    log_folder = ptt.join(daily_dir, "logs", "Summary")
    if control:
        log_folder = ptt.join(log_folder, 'control')
    if setup.epoch:
        log_folder = ptt.join(log_folder, 'epoch')
        setup.epoch = True
    if setup.custom is not None:
        log_folder = ptt.join(log_folder,setup.custom)
    os.makedirs(log_folder, exist_ok = True)
    return(log_folder)

def cleanup_bkups(setup):
    allsky = False if setup.custom is not None else True
    fnames.set(setup.boss_spectro_redux, setup.run2d, dev=False,
                           epoch=setup.epoch, allsky=allsky, custom=setup.custom)
    bk_files = OrderedDict()
    fnames.bk.set(flag = '*')
    for bf in glob(fnames.bk.spAllfile):
        bk_files[bf] = bf.split('-')[-1]
    idx = np.argsort(np.asarray(list(bk_files.values())))
    for i, id in enumerate(np.flip(idx)):
        key = list(bk_files.keys())[id]
        if i > setup.backup -1:
            splog.debug(f'Removing Backup: {key} ({i+1})')
            os.remove(key)
            splog.debug(' '*16+key.replace('spAll-','spAll-lite-'))
            os.remove(key.replace('spAll-','spAll-lite-'))
            splog.debug(' '*16+key.replace('spAll-','spAllLine-'))
            os.remove(key.replace('spAll-','spAllLine-'))
        else:
            splog.debug(f'Keeping Backup: {key} ({i+1})')
    return

def build(setup, no_submit=False, daily=False,
            email_start = False, obs = None):
    if not noslurm:
        queue1 = queue()
        queue1.verbose = True
    else:
        queue1 = None
    log_folder = _build_log_dir(setup, control = True)
    os.makedirs(ptt.join(log_folder), exist_ok = True)

    log = ptt.join(_build_log_dir(setup, control = False), setup.run2d, "pySummary_"+jdate.astype(str))

    if not daily:
        splog.open(ptt.join(log_folder, jdate.astype(str)+'.log'))
    else:
        splog.add_file(ptt.join(log_folder, jdate.astype(str)+'.log'))
    
    if email_start:
        splog.emailer()

    
    if setup.daily is True:
        with fits.open(ptt.join(setup.boss_spectro_redux, setup.run2d,
                                'spAll-'+setup.run2d+'.fits')) as ff:
            tf = Table(ff[1].data)
            latest_mjd = tf['MJD'].max()

        if not check_daily(setup.module, os.getenv('DAILY_DIR'), latest_mjd):
            splog.debug('Skipping run')
            splog.send_email('fieldmerge '+args.run2d +' MJD='+jdate.astype(str),
                      ptt.join(os.getenv('DAILY_DIR'), 'etc','emails'))
            return()
        if not check_fieldlist(setup.boss_spectro_redux, setup.run2d, latest_mjd):
            splog.debug('SpAll-'+setup.run2d+' up to date')
            splog.send_email('fieldmerge '+setup.run2d +' MJD='+jdate.astype(str),
                      ptt.join(os.getenv('DAILY_DIR'), 'etc','emails'))
            return()
    
    
    splog.info('===============================================')
    splog.debug(time.ctime())

    splog.info(setup)

    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout



    title = setup.run2d+'/apo_lco/'+jdate.astype(str)+'/BOSS_Summary'
    if setup.epoch: title = title+'/epoch'
    if setup.custom is not None:
        title = title+'/'+setup.custom
    
    if setup.nodes > 1:
        title = title.replace('/','_')


    if not noslurm:
        queue1.create(label = title, nodes = setup.nodes, ppn = setup.ppn,
                      walltime = setup.walltime,
                      alloc = setup.alloc, partition = setup.partition,
                      mem_per_cpu = setup.mem_per_cpu, shared = setup.shared)


    job_dir = ptt.join(os.getenv('SLURM_SCRATCH_DIR'),title,queue1.key)
    full_cmd = ["#!/usr/bin/env bash"]
    full_cmd.append(f"cd {ptt.abspath(ptt.join(log_folder,'..'))}")
    if setup.module is not None:
        full_cmd.append("module purge ; module load "+setup.module)
    full_cmd.append("set -o verbose")
    if (setup.custom is None):
        if not setup.no_fieldlist:
            full_cmd.append(f"fieldlist --create --run1d {setup.run2d} --run2d {setup.run2d}")
    fm_cmd = f"fieldmerge --lite --include_bad --XCSAO"
    bk_cmd = f"cleanup_backups --topdir {setup.boss_spectro_redux} --run2d {setup.run2d}"
    if setup.merge_only:
        fm_cmd = fm_cmd+" --merge_only"
    if setup.limit is not None:
        fm_cmd = fm_cmd+f" --limit {setup.limit}"
    if setup.backup is not None:
        fm_cmd = fm_cmd+" --bkup"
        bk_cmd = bk_cmd+f" --backups {setup.backup}"
    if setup.epoch:
        fm_cmd = fm_cmd+" --epoch"
        bk_cmd = bk_cmd+" --epoch"
    if setup.custom is not None:
        bk_cmd = bk_cmd+" --custom {setup.custom}"
        fm_cmd = fm_cmd+f" --allsky --custom {setup.custom}"
    if setup.skip_specprimary:
        fm_cmd = fm_cmd+" --skip_specprimary"
    if setup.update_specprimary:
        fm_cmd = fm_cmd+" --update_specprimary"
    if setup.verbose:
        fm_cmd = fm_cmd+" --verbose"
    if setup.n_iter is None:
        setup.n_iter = 1
    if setup.ndays is not None:
        fm_cmd = fm_cmd+f" --ndays {setup.ndays}"
    
    dlog_folder = _build_log_dir(setup, control = False)

    
    for i in range(setup.n_iter):
        if setup.backup is not None:
            if i > 0:
                full_cmd.append(f'current_time=$(date "+%Y.%m.%d-%H.%M.%S")')
                full_cmd.append(f"cp -p {log}.o.log {log}.o.log-$current_time ")
                full_cmd.append(f"cp -p {log}.e.log {log}.o.log-$current_time ")
            full_cmd.append("")
        if setup.n_iter == 1:
            full_cmd.append(fm_cmd)
        else:
            full_cmd.append(fm_cmd)
        
    if setup.backup is not None:
        full_cmd.append(bk_cmd)
        
    if setup.utah_daily:
        full_cmd.append(f'sas_mos_too boss -t apo25m -d {setup.run2d}')
        full_cmd.append(f'sas_mos_too boss -t lco25m -d {setup.run2d}')
        full_cmd.append(f'sdss5db_update_boss -d {setup.run2d} -v -s -Y -p')
        
    with open(ptt.join(job_dir,'run_pySummary'),'w') as r:
        for c in full_cmd:
            r.write(c+'\n')
    if not noslurm:
        queue1.append("source "+ptt.join(job_dir,'run_pySummary'),
                      outfile = log+".o.log", errfile = log+".e.log")
    else:
        print("source "+ptt.join(job_dir,'run_pySummary')+'>'+log+".o.log 2> "+log+".e.log")
    if obs is not None:
        obs = np.atleast_1d(obs)
        lcoflag = ' --lco' if obs[0].upper() == 'LCO' else ''
        epochflag = ' --epoch' if setup.epoch else ''
        if not noslurm:
            queue1.append(f"plot_QA    --run2d {setup.run2d} {lcoflag} {epochflag} ; ")
        else:
            print(f"plot_QA    --run2d {setup.run2d} {lcoflag} {epochflag} ; ")
    

    if not noslurm:
        queue1.commit(submit=(not no_submit))

    output = new_stdout.getvalue()
    sys.stdout = old_stdout
    splog.info(output)

    
    subject = _build_subject(setup, jdate.astype(str))
    
    if email_start:
        splog.send_email(subject, ptt.join(os.getenv('DAILY_DIR'), 'etc','emails'))
    return(queue1, title, [log+".o.log",log+".e.log"])


