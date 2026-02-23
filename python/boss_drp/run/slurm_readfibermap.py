#!/usr/bin/env python3
from boss_drp.utils import load_env
from boss_drp.field import field_to_string
from boss_drp.field import Field
from boss_drp import daily_dir
from boss_drp.Config import config
from boss_drp.run.queue import Queue
import sys

import argparse
from os import getenv, makedirs
import os.path as ptt
from pydl.pydlutils.yanny import read_table_yanny
from glob import glob


def setup_run():
    config.add_config('readfibermap')
       

def slurm_readfibermap():

    setup_run()
    
    if ppn is not None:
        maxppn = min(20,int(config.readfibermap_queue.get('ppn')))
        if ppn > int(config.readfibermap_queue.get('ppn')):
            print(f'Resetting to {maxppn} due to cluster configuration')
            ppn = int(config.readfibermap_queue.get('ppn'))
        if ppn > 20:
            if input('Are you sure you want to run more then 20 in parallel? (y/[n])  ').lower() == 'y':
                pass
            else:
                ppn = maxppn
                print(f'Reseting to {maxppn} in parallel')
        config.readfibermap_queue.set('ppn', ppn)
    
    
    fdir = Field(config.pipe['general.BOSS_SPECTRO_REDUX'],config.pipe['general.RUN2D'], '*').dir()
    plan2ds = glob(ptt.join(fdir,'spPlan2d*.par'))
    
    qu = build(plan2ds, daily_dir=daily_dir,
                mjd=config.pipe['fmjdselect.mjd'],
                mjdstart= config.pipe['fmjdselect.mjdstart'], 
                mjdend=config.pipe['fmjdselect.mjdend'], 
                obs = config.pipe['fmjdselect.obs'])

def build(plan2ds, daily=False, 
            mjd=None, mjdstart= None, mjdend=None, no_submit=False,
            obs = ['apo','lco'], daily_dir=daily_dir):
    i = 0
    clobber = config.pipe['Clobber.clobber_fibermap']
    V_TARG = config.pipe['general.V_TARG']
    title = 'readfibermap_'+config.pipe['general.RUN2D']
    #if not daily:
    log = ptt.join(daily_dir, "logs", "readfibermap", config.pipe['general.RUN2D'], "readfibermap_")
    makedirs(ptt.join(daily_dir, "logs", "readfibermap", config.pipe['general.RUN2D']), exist_ok = True)
    cmds = []
    for plan2d in plan2ds:
        thisplan = read_table_yanny(plan2d, 'SPEXP')
        thisplan.convert_bytestring_to_unicode()
        thismjd = int(thisplan.meta['MJD'])
        if mjd is not None:
            if thismjd not in mjd:
                continue
        else:
            if mjdstart is not None:
                if thismjd < mjdstart:
                    continue
            if mjdend is not None:
                if thismjd > mjdend:
                    continue
        try:
            if thisplan.meta['OBS'].lower() not in obs:
                continue
        except:
            if thisplan['name'][0][0].split('-')[1] in ['b2','r2']:
                tobs = 'LCO'
            else:
                tobs = 'APO'
            if tobs.lower() not in obs:
                continue
        if not clobber:
            if ptt.exists(plan2d.replace('spPlan2d','spfibermap').replace('.par','.fits')):
                continue
        
        thislog = log+ptt.basename(plan2d).replace('spPlan2d-','').replace('.par','')
        drf = '' if V_TARG == '*' else f' --V_TARG {V_TARG}'
        thiscmd = (f"cd {ptt.dirname(plan2d)} ; " +
                  f"readfibermaps --spplan2d {ptt.basename(plan2d)} {drf}")
                
        if clobber:
            thiscmd = thiscmd+' --clobber'
        cmds.append((thiscmd,thislog))

    if len(cmds) > 0:
        if len(cmds) < config.readfibermap_queue.get('ppn'):
            config.readfibermap_queue.set('ppn', max([len(cmds), 2]))

        print(config.readfibermap_queue.to_str())
        queue1 = Queue(config.readfibermap_queue, key=None, verbose=True)
        queue1.create(**config.readfibermap_queue.to_dict(label=title))


        for (thiscmd, thislog) in cmds:
            queue1.append(thiscmd, outfile = thislog+".o.log", errfile = thislog+".e.log")
 
        queue1.commit(hard=True,submit=not config.readfibermap_queue.get('no_submit'))
    
        return(queue1)
    else:
        print('No fibermaps Built')
    return(None)

