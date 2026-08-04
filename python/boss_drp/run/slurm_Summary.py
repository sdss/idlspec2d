#!/usr/bin/env python3
from jinja2 import Template
from boss_drp.post.fieldmerge import fieldlist_name
from boss_drp.summary import Summary_names
from boss_drp.utils import jdate, send_email
from boss_drp import daily_dir
from boss_drp.utils.splog import splog
from boss_drp.Config import config, fill_none_with_false
from boss_drp.run.queue import Queue
import boss_drp
from pydl.pydlutils.yanny import yanny
from astropy.io import fits
from astropy.table import Table

import os
import os.path as ptt
import numpy as np
import time
import re


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
    if ptt.exists(fieldlist_name.parquet):
        flist = Table.read(fieldlist_name.parquet)
    elif ptt.exists(fieldlist_name.name):
        flist = Table.read(fieldlist_name.name)
    else:
        return(False)
    r = re.compile('Done[\w]*', re.IGNORECASE)
    idx = [i for i, x in enumerate(flist['STATUS1D'].data ) if r.search(x)]
    flist = flist[idx]
    latest_redux = max(flist['MJD'])
    return(latest_redux > spall_mjd)


def Setup():
    config.add_config('Summary')
    fill_none_with_false(config.Summary_queue)
def slurm_Summary():
    Setup()
 
    queue1, title, attachements = build()#setup, email_start = email_start)
    monitor = config.pipe['monitor.pipe_monitor']
    monitor  = queue1.monitor_job(monitor=monitor, pause = 300, jobname = title)                   
                          
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
                        if 'Successful completion of build_spall' in line:
                            flags.append('Complete: fieldmerge')
                        elif 'Successful completion of fieldmerge' in line:
                            flags.append('Complete: fieldmerge')
                        elif 'Successful completion of fieldlist' in line:
                            flags.append('Complete: fieldlist')
                        elif 'exception:' in line:
                            flags.append('Crashed')
                        elif 'EXITING' in line:
                            flags.append('Killed')
                        elif 'Killed' in line:
                            flags.append('Killed')
                        elif 'Traceback' in line:
                            flags.append('Crashed')
                        elif 'errno' in line:
                            flags.append('Crashed')
            if 'Killed' in flags:
                subject = 'Killed: '+subject
            elif 'Crashed' in flags:
                subject = 'Crashed: '+subject
            elif ('Complete: fieldmerge' in flags and 'Complete: fieldlist' in flags):
                subject = 'Complete: '+subject
            elif 'Complete: fieldmerge' in flags and (not config.pipe['Stage.run_fieldlist']):
                subject = 'Complete: '+subject
            elif 'Complete: fieldmerge' in flags and 'Complete: fieldlist' not in flags:
                subject = 'Incomplete: '+subject +' (fieldlist not complete)'
            elif 'Complete: fieldlist' in flags and 'Complete: fieldmerge' not in flags:
                subject = 'Incomplete: '+subject +' (fieldmerge not complete)'
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
    elif config.pipe['customSettings.custom_name'] is not None:
        log_folder = ptt.join(log_folder,config.pipe['customSettings.custom_name'])
    else:
        log_folder = ptt.join(log_folder, 'daily')
    os.makedirs(log_folder, exist_ok = True)
    return(log_folder)

def build():#setup, daily=False, email_start = False, obs = None):
    log_folder = _build_log_dir(control = True)
    dlog_folder = _build_log_dir(control = False)
    obs = config.pipe['fmjdselect.obs'] or None
    if isinstance(obs, list):
        if len(obs) == 2:
            obs = None

    os.makedirs(ptt.join(log_folder), exist_ok = True)

    log = ptt.join(_build_log_dir(control = False), config.pipe['general.RUN2D'], "pySummary_"+jdate.astype(str))

    if not config.pipe['Sumamry.batchwise.after_daily']:
        splog.open(ptt.join(log_folder, jdate.astype(str)+'.log'))
    else:
        splog.add_file(ptt.join(log_folder, jdate.astype(str)+'.log'))
    
    if config.pipe['email.email_start']:
        splog.emailer()

    
    if config.pipe['Sumamry.batchwise.after_daily']:
        summ = Summary_names()
        summ.set(indir = config.pipe['general.BOSS_SPECTRO_REDUX'], 
                 run2d = config.pipe['general.RUN2D'])
        if ptt.exists(summ.spAllfile_parquet):
            spall = Table.read(summ.spAllfile_parquet)
            latest_mjd = spall['MJD'].max()
            spall = None            
        elif ptt.exists(summ.spAllfile):
            spall = Table.read(summ.spAllfile)
            latest_mjd = spall['MJD'].max()
            spall = None
        else:
            splog.debug('No spAll file found')
            latest_mjd = 0

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

    if obs is not None:
        obsstr = '_'.join(np.atleast_1d(obs).tolist()).upper()
    else:
        obsstr = 'apo_lco'
    title = config.pipe['general.RUN2D']+f'/{obsstr}/'+jdate.astype(str)+'/BOSS_Summary'
    if config.pipe['fmjdselect.epoch']: title = title + '/epoch'
    if config.pipe['customSettings.custom_name'] is not None:
        title = title+'/'+config.pipe['customSettings.custom_name']
    
    if config.Summary_queue.get('nodes') > 1:
        title = title.replace('/','_')

    with splog.capture_prints():

        queue1 = Queue(config.Summary_queue, verbose=True)
        queue1.create(**config.Summary_queue.to_dict(label=title))
        job_dir = ptt.join(config.Summary_queue.get('queue_sub_dir',os.getcwd()),title,queue1.key)

        flags = []
        bkflags = []
        for key, value in config.pipe['Summary'].items():
            if key == 'fieldlist':
                continue
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
            if value is None:
                continue
            flags.append(f'--{key} {value}')

        key_map = {'backup':'bkup', 'clobber_fmjd':'clobber'}
        for key, value in config.pipe['Summary.batchwise'].items():
            if key in ['database']:
                continue
            keym = key_map[key] if key in key_map else key

            if isinstance(value, bool) or str(value).lower() in ['true', 'false']:
                if str(value).lower() == 'true':
                    flags.append(f'--{keym}')
                continue
            if keym in ['bkup']:
                flags.append(f'--bkup')
                continue
            if isinstance(value, dict):
                continue
            if value is None: 
                continue
            flags.append(f'--{keym} {value}')

        if config.pipe['fmjdselect.epoch']:
            flags.append('--epoch')
            bkflags.append('--epoch')
        if config.pipe['customSettings.custom_name'] is not None:
            flags.append(f"--custom {config.pipe['customSettings.custom_name']}")
            bkflags.append(f"--custom {config.pipe['customSettings.custom_name']}")
        if config.pipe['customSettings.allsky']:
            flags.append(f"--allsky")
        if config.pipe['Summary.batchwise.backup'] is not None:
            bkflags.append(f"--backups {config.pipe['Summary.batchwise.backup']}")

        if config.pipe['general.BOSS_SPECTRO_REDUX'] != os.getenv('BOSS_SPECTRO_REDUX'):
            flags.append(f"--indir {config.pipe['general.BOSS_SPECTRO_REDUX']}")
            bkflags.append(f"--topdir {config.pipe['general.BOSS_SPECTRO_REDUX']}")

        if config.pipe['general.RUN2D'] != os.getenv('RUN2D'):
            flags.append(f"--run2d {config.pipe['general.RUN2D']}")
            bkflags.append(f"--run2d {config.pipe['general.RUN2D']}")

        if config.pipe['general.RUN1D'] != os.getenv('RUN1D'):
            flags.append(f"--run1d {config.pipe['general.RUN1D']}")

        fieldmergeflags = ' '.join(flags)
        bkflags = ' '.join(bkflags)

        fieldmergeflags_itter = fieldmergeflags.replace(' --lite','')
        bk_cmd = (f"boss_drp clean backups {bkflags}")

        flist_flags = []
        if config.pipe['Stage.run_fieldlist']:
            flist_flags.append('--create')
            if config.pipe['general.BOSS_SPECTRO_REDUX'] != os.getenv('BOSS_SPECTRO_REDUX'):
                flist_flags.append('--topdir '+config.pipe['general.BOSS_SPECTRO_REDUX'])
            flist_flags.append('--run1d '+config.pipe['general.RUN1D'])
            flist_flags.append('--run2d '+config.pipe['general.RUN2D'])
            if config.pipe['fmjdselect.epoch']:
                flist_flags.append('--epoch')
            for key, value in config.pipe['Summary.fieldlist'].items():
                if isinstance(value, bool) or str(value).lower() in ['true', 'false']:
                    if str(value).lower() == 'true':
                        flags.append(f'--{key}')
                    continue
                if isinstance(value, dict):
                    continue
                if value is None:
                    continue
                flags.append(f'--{key} {value}')
                
        pipe_flags = dict(control_dir = ptt.abspath(ptt.join(log_folder,'..')),
                        module = config.pipe['general.module'], RUN2D=config.pipe['general.RUN2D'],
                        log = log, flist_flags = ' '.join(flist_flags),
                        fieldlist = config.pipe['Stage.run_fieldlist'],
                        n_iter = config.pipe['Summary.batchwise.n_iter'] or 1,
                        fieldmergeflags = fieldmergeflags,
                        fieldmergeflags_itter=fieldmergeflags_itter,
                        bk_cmd = bk_cmd,
                        database = config.pipe['Summary.batchwise.database'] 
                        )


        template = ptt.join(ptt.dirname(boss_drp.__file__), 'etc','templates','Summary.j2')
        with open(ptt.join(job_dir,'run_pySummary'), "w", encoding="utf-8") as output_file:
            with open(template) as template_file:
                j2_template = Template(template_file.read())
                output_file.write(re.sub(r'\n\s*\n+', '\n\n',j2_template.render(pipe_flags)))
        splog.info("Generated Summary script: %s", ptt.join(job_dir,'run_pySummary'))
        queue1.append("source "+ptt.join(job_dir,'run_pySummary'),
                        outfile = log+".o.log", errfile = log+".e.log")
        if obs is not None:
            obs = np.atleast_1d(obs)
            lcoflag = ' --lco' if obs[0].upper() == 'LCO' else ''
            epochflag = ' --epoch' if config.pipe['fmjdselect.epoch'] else ''
            queue1.append(f"plot_QA    --run2d {config.pipe['general.RUN2D']} {lcoflag} {epochflag} ; ")
        
        queue1.commit(submit=(not config.Summary_queue.get('no_submit')))
    
    subject = _build_subject(jdate.astype(str))
    
    if config.pipe['email.email_start']:
        splog.send_email(subject, ptt.join(daily_dir, 'etc','emails'))
    return(queue1, title, [log+".o.log",log+".e.log"])


#TODO: MJD range???