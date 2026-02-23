#!/usr/bin/env python3
from boss_drp.utils import get_dirs, load_env, mjd_match
from boss_drp.field import Field, field_to_string
from boss_drp.run.config2redux import config2redux
from boss_drp.utils import jdate
from boss_drp import daily_dir
from boss_drp.utils.splog import splog
from boss_drp.Config import config
from boss_drp.run.queue import Queue
import argparse
from os import getenv
import os.path as ptt
from glob import glob
from pydl.pydlutils.yanny import yanny, read_table_yanny
import sys
import numpy as np
import io
import datetime
import astropy.time

def uubatchpbs( daily=False):


    if daily is False:
        splog.emailer()
        
    error = False
    
    if not ptt.isdir(config.pipe['general.BOSS_SPECTRO_REDUX']):
        splog.warning('topdir (BOSS_SPECTRO_REDUX) is invalid')
        error = True
    else:
        splog.info('topdir   '+config.pipe['general.BOSS_SPECTRO_REDUX'])

    splog.info('Config:\n'+str(config))
    
    obs = config.pipe['fmjdselect.obs']

    afc = Field(config.pipe['general.BOSS_SPECTRO_REDUX'],
                config.pipe['general.RUN2D'], '*')
    if not config.pipe['customSettings.allsky']:
        fielddirs = get_dirs(ptt.dirname(afc.dir()), field=True,
                             match = config.pipe['fmjdselect.field'],
                             start = config.pipe['fmjdselect.fieldstart'], 
                             end   = config.pipe['fmjdselect.fieldend'])
        if len(fielddirs) == 0:
            splog.warning('No Directories Found')
            error = True
    else:
        if len(obs) == 2:
            fielddirs = [config.pipe['customSettings.custom_name']]
        else:
            fielddirs = []
        fielddirs.extend([config.pipe['customSettings.custom_name']+'_'+ob for ob in obs])
        
    if error:
        splog.close_elogger()
        return()

    redux_list = []
    skipped = 0
   
    if config.pipe['fmjdselect.epoch'] is False:
        if config.pipe['customSettings.custom_name'] is None:
            plan_str = 'spPlancomb-*.par'
        else:
            plan_str = 'spPlanCustom-{custom}-*.par'
            plan_str_bkup = plan_str
    else:
        if custom is None:
            plan_str = 'spPlancombepoch-*.par'
        else:
            plan_str = 'spPlancombepoch_{custom}-*.par'
            plan_str_bkup = plan_str

    if config.pipe['Clobber.clobber_pipe']:
        if ((config.pipe['fmjdselect.mjd'] is None) and 
            (config.pipe['fmjdselect.mjdstart'] is None) and 
            (config.pipe['fmjdselect.mjdend'] is None)):
            splog.info('No MJDs Selected while clobber is set')
            val = input('Do you want to continue? (yes/NO)')
            if val.lower() != 'yes':
                exit()

    for fielddir in fielddirs:
        if config.pipe['customSettings.custom_name'] is not None:
            plan_str = plan_str_bkup
            plan_str = plan_str.format(custom=fielddir)
        fc = Field(config.pipe['general.BOSS_SPECTRO_REDUX'], 
                   config.pipe['general.RUN2D'], fielddir, 
                   custom_name = config.pipe['customSettings.custom_name'], 
                   epoch = config.pipe['fmjdselect.epoch'])
        planfile = glob(ptt.join(fc.dir(), plan_str))
        for plancombine in planfile:
            yplan = yanny(plancombine)
            hdr = yplan.new_dict_from_pairs()
            thismjd = hdr['MJD']
            if 'OBS' in hdr.keys():
                thisobs = hdr['OBS'].lower()
                if thisobs.lower() not in obs:
                    continue
            else:
                thisobs = None
            if not mjd_match(thismjd, mjd=config.pipe['fmjdselect.mjd'], 
                                mjdstart=config.pipe['fmjdselect.mjdstart'], 
                                mjdend=  config.pipe['fmjdselect.mjdend']):
                continue    
            if config.pipe['customSettings.custom_name'] is None:
                plan2d = hdr['planfile2d'].replace("'",'').split(' ')
                field = hdr['fieldid']
                mjd = hdr['MJD']
            else:
                plan2d = []
                if config.pipe['customSettings.allsky']:
                    field = None
                    mjd = hdr['CreateMJD']
                else:
                    field = hdr['fieldid']
                    mjd = hdr['CreateMJD']
            
            cmjds = [None]
            if config.pipe['customSettings.custom_single_mjd']:
                cmjds = yplan['COADDPLAN']['EPOCH_COMBINE'] 

            for cmjd in set(cmjds):
                redux1 = config2redux(plan2d, plancombine, obs=thisobs,
                                    field = field, mjd = mjd, 
                                    custom = config.pipe['customSettings.custom_name'],
                                    allsky=config.pipe['customSettings.allsky'], 
                                    custom_single_mjd=cmjd)
                if redux1 is None:
                    skipped += 1
                    continue
                redux_list.append(redux1)

    splog.info('')
    splog.info('---------------------------------------------------')
    splog.info('boss_redux: #BOSS Field-MJDs Done  = '+str(skipped))
    splog.info('            #BOSS Field-MJDs To Do = '+str(len(redux_list)))
    splog.info('---------------------------------------------------')
    splog.info('')
    
    if len(redux_list) == 0: 
        splog.close_elogger()
        return None, redux_list
    
    if config.queue.get('max.ppn') is not None:
        max_c = int(config.queue.get('max.ppn'))
    else:
        max_c = 1

    if config.queue.get('ppn') is None:
        config.queue.set('ppn', np.min([max_c, np.max([len(redux_list),2])]))
    
    #Checking fallback condition - some queues have strict requirements 
    config.queue.check_fallback()

    if daily is True:
        if config.queue.get('nodes') is None and len(redux_list) > config.queue.get('ppn'):
            config.queue.set('nodes', int(np.ceil(len(redux_list)/max_c)))
        else:
            config.queue.set('nodes', 1)
    elif config.queue.get('nodes') is None:
        config.queue.set('nodes',1)
        
    run2d = config.pipe['general.RUN2D']
    custom = config.pipe['customSettings.custom_name']
    epoch = config.pipe['fmjdselect.epoch']

    obsstr = '_'.join(obs).upper()
    if daily: 
        mjdstr = str(mjd[0]) 
    elif len(mjd) == 1:
        mjdstr = str(mjd[0])
    else:
        mjdstr = 'batch'
    if custom is None:
        if not epoch:
            label = f'{run2d}/{obsstr}/{mjdstr}/daily/'
        else:
            label = f'{run2d}/{obsstr}/{mjdstr}/epoch/'
    else:
        if not epoch:
            label = f'{run2d}/{obsstr}/{mjdstr}/{custom}/'
        else:
            label = f'{run2d}/{obsstr}/{mjdstr}/epoch/{custom}/'


    if config.queue.get('nodes') > 1:
        label = label.replace('/','_')

    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    queue1 = Queue(config.queue, key=None, verbose=True)
    queue1.create(**config.queue.to_dict(label=label))
    rlist = redux_list.copy()
    for i in range(len(redux_list)):
        cmd, log, err = make_run_cmd(redux_list[0])
        redux_list.pop(0)
        queue1.append(cmd, outfile = log, errfile = err)
    queue1.commit(hard=True, submit=(not config.queue.get('no_submit')))
    output = new_stdout.getvalue()
    sys.stdout = old_stdout
    splog.info(output)

    if config.pipe['email.sendemail'] is True:
        if daily:
            tjdate = mjd[0]
        else:
            tjdate = jdate.astype(str)
        splog.send_email('UUBATCH '+run2d +' MJD='+str(tjdate) +' OBS='+','.join(obs),
                              ptt.join(daily_dir, 'etc','emails'), 
                              allemail=config.pipe['email.allemail'])
    else:
        splog.close_elogger()

    return(queue1, rlist)

def make_run_cmd(redux):
    cmd = 'source '+redux
    log = redux+'.o'
    err = redux+'.e'
    return(cmd, log, err)

