#!/usr/bin/env python3
from jinja2 import Template
from boss_drp.post.fieldmerge import summary_names as fnames, fieldlist_name
from boss_drp.utils import jdate, send_email
from boss_drp import daily_dir
from boss_drp.utils.splog import splog
from boss_drp.Config import config
from boss_drp.run.queue import Queue
import boss_drp
from pydl.pydlutils.yanny import yanny
from astropy.io import fits
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


def Setup():
    config.add_config('Summary')

def slurm_Summary():
    Setup()
 
    
<<<<<<< Updated upstream
    if full:
        setup.shared = False
        setup.ppn = os.getenv('SLURM_PPN')
        setup.mem = 0 #500000
        setup.mem_per_cpu = None
    else:
        setup.ppn = 10
        if mem is not None:
            setup.mem = mem
            setup.mem_per_cpu = None
            #setup.mem_per_cpu  = mem/setup.ppn
    setup.walltime = walltime
    
    queue1, title, attachements = build(setup, no_submit=no_submit,
                                        email_start = email_start)
                                
                                
    if no_submit:
        monitor = False
=======
    queue1, title, attachements = build()#setup, email_start = email_start)
    monitor = config.pipe['monitor.pipe_monitor']
    monitor  = queue1.monitor_job(monitor=monitor, pause = 300, jobname = title)                   
                          
>>>>>>> Stashed changes
    if monitor:
                
        subject = _build_subject(jdate.astype(str))
        splog.close_file()
        flags = []
        if config.pipe['Stage.run_fieldlist']:
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
                      attachements)

def _build_subject(mjd):
    mstr = config.pipe['general.module'] if config.pipe['general.module'] is not None else config.pipe['general.RUN2D']
    if config.pipe['fmjdselect.epoch']:
        subject = f'BOSS Summary {mstr} epoch MJD={mjd}'
    elif config.pipe['customSettings.customname'] is not None:
        subject = f"BOSS Summary {mstr} {config.pipe['customSettings.customname']} MJD={mjd}"
    else:
        subject = f'BOSS Summary {mstr}  MJD={mjd}'
    return(subject)

def _build_log_dir(control = False):
    log_folder = ptt.join(daily_dir, "logs", "Summary")
    if control:
        log_folder = ptt.join(log_folder, 'control')
    if config.pipe['fmjdselect.epoch']:
        log_folder = ptt.join(log_folder, 'epoch')
    elif config.pipe['customSettings.customname'] is not None:
        log_folder = ptt.join(log_folder,config.pipe['customSettings.customname'])
    else:
        log_folder = ptt.join(log_folder, 'daily')
    os.makedirs(log_folder, exist_ok = True)
    return(log_folder)

def build():#setup, daily=False, email_start = False, obs = None):
    log_folder = _build_log_dir(control = True)
    dlog_folder = _build_log_dir(control = False)
    obs = config.pipe['fmjdselect.obs']

    os.makedirs(ptt.join(log_folder), exist_ok = True)

    log = ptt.join(_build_log_dir(control = False), config.pipe['general.RUN2D'], "pySummary_"+jdate.astype(str))

    if not config.pipe['Sumamry.batchwise.after_daily']:
        splog.open(ptt.join(log_folder, jdate.astype(str)+'.log'))
    else:
        splog.add_file(ptt.join(log_folder, jdate.astype(str)+'.log'))
    
    if config.pipe['email.email_start']:
        splog.emailer()

    
    if config.pipe['Sumamry.batchwise.after_daily']:
        with fits.open(ptt.join(config.pipe['general.BOSS_SPECTRO_REDUX'], 
                                config.pipe['general.RUN2D'],
                                'spAll-'+config.pipe['general.RUN2D']+'.fits')) as ff:
            tf = Table(ff[1].data)
            latest_mjd = tf['MJD'].max()

        if not check_daily(config.pipe['general.module'], daily_dir, latest_mjd):
            splog.debug('Skipping run')
            splog.send_email('fieldmerge '+config.pipe['general.RUN2D'] +' MJD='+jdate.astype(str),
                      ptt.join(daily_dir, 'etc','emails'))
            return()
        if not check_fieldlist(config.pipe['general.BOSS_SPECTRO_REDUX'],
                                config.pipe['general.RUN2D'], latest_mjd):
            splog.debug('SpAll-'+config.pipe['general.RUN2D']+' up to date')
            splog.send_email('fieldmerge '+config.pipe['general.RUN2D'] +' MJD='+jdate.astype(str),
                      ptt.join(daily_dir, 'etc','emails'))
            return()
    
    
    splog.info('===============================================')
    splog.debug(time.ctime())

    splog.info(config)

    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout



    title = config.pipe['general.RUN2D']+'/apo_lco/'+jdate.astype(str)+'/BOSS_Summary'
    if config.pipe['fmjdselect.epoch']: title = title+'/epoch'
    if config.pipe['customSettings.custom_name'] is not None:
        title = title+'/'+config.pipe['customSettings.custom_name']
    
    if config.Summary_queue.get('nodes') > 1:
        title = title.replace('/','_')

    queue1 = Queue(config.Summary_queue, verbose=True)

    queue1.create(**config.Summary_queue.to_dict(label=title))
    job_dir = ptt.join(config.Summary_queue.get('queue_sub_dir',os.getcwd()),title,queue1.key)



    flags = []
    bkflags = []
    for key, value in config.pipe['Summary'].items():
        if isinstance(value, bool) or str(value).lower() in ['true', 'false']:
            if str(value).lower() == 'true':
                flags.append(f'--{key}')
            continue
        if isinstance(value, dict):
            continue
        if key == 'skip_specprimary':
            if value == 'update':
                flags.append('--update_specprimary')
                continue
        flags.append(f'--{key} {value}')

    key_map = {'backup':'bkup'}
    for key, value in config.pipe['Summary.batchwise'].items():
        if key in ['database']:
            continue
        keym = key_map[key] if key in key_map else key

        if isinstance(value, bool) or str(value).lower() in ['true', 'false']:
            if str(value).lower() == 'true':
                flags.append(f'--{keym}')
            continue
        if isinstance(value, dict):
            continue
        flags.append(f'--{keym} {value}')

    if config.pipe['fmjdselect.epoch']:
        flags.append('--epoch')
        bkflags.append('--epoch')
    if config.pipe['customSettings.customname']:
        flags.append(f"--custom {config.pipe['customSettings.customname']}")
        bkflags.append(f"--custom {config.pipe['customSettings.customname']}")
    if config.pipe['customSettings.allsky']:
        flags.append(f"--allsky")
    if config.pipe['Summary.batchwise.backup']:
        bkflags.append(f"--backups {config.pipe['Summary.batchwise.backup']}")

    fieldmergeflags = ' '.join(flags)
    bkflags = ' '.join(bkflags)

    fieldmergeflags_itter = fieldmergeflags.replace(' --lite','')
    bk_cmd = (f"cleanup_backups --topdir {config.pipe['general.BOSS_SPECTRO_REDUX']} "+
              f"--run2d {config.pipe['general.RUN2D']} {bkflags}")


            
    pipe_flags = dict(control_dir = ptt.abspath(ptt.join(log_folder,'..')),
                      module = config.pipe['general.module'], RUN2D=config.pipe['general.RUN2D'],
                      RUN1D = config.pipe['general.RUN1D'], log = log,
                      fieldlist = config.pipe['stage.run_fieldlist'],
                      n_iter = config.pipe['Summary.batchwise.n_iter'],
                      fieldmergeflags = fieldmergeflags,
                      fieldmergeflags_itter=fieldmergeflags_itter,
                      bk_cmd = bk_cmd,
                      database = config.pipe['Summary.batchwise.database'] 
                     )


    template = ptt.join(ptt.dirname(boss_drp.__file__), 'etc','templates','Summary.j2')
    with open(ptt.join(job_dir,'run_pySummary'), "w", encoding="utf-8") as output_file:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            output_file.write(j2_template.render(pipe_flags))
    
    queue1.append("source "+ptt.join(job_dir,'run_pySummary'),
                      outfile = log+".o.log", errfile = log+".e.log")
    if obs is not None:
        obs = np.atleast_1d(obs)
        lcoflag = ' --lco' if obs[0].upper() == 'LCO' else ''
        epochflag = ' --epoch' if config.pipe['fmjdselect.epoch'] else ''
        queue1.append(f"plot_QA    --run2d {config.pipe['general.RUN2D']} {lcoflag} {epochflag} ; ")
    
    queue1.commit(submit=(not config.Summary_queue.get('no_submit')))

    output = new_stdout.getvalue()
    sys.stdout = old_stdout
    splog.info(output)

    
    subject = _build_subject(jdate.astype(str))
    
    if config.pipe['email.email_start']:
        splog.send_email(subject, ptt.join(daily_dir, 'etc','emails'))
    return(queue1, title, [log+".o.log",log+".e.log"])


