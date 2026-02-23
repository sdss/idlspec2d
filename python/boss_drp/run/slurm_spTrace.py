#!/usr/bin/env python3
from boss_drp.prep.spplan_trace import spplanTrace
from boss_drp.utils import load_env
from boss_drp.Config import config
from boss_drp.run.queue import Queue
import sys
from os import path as ptt
import numpy as np
import datetime


run2d = None
topdir = None

def setup_run(nodes=None, alloc=None, partition=None, nbundle=None, 
              mjd = [], maxjobs=None):
    config.add_config('spTrace')

    config.spTrace_queue.set('nodes',nodes)
    config.spTrace_queue.set('alloc',alloc)
    config.spTrace_queue.set('partition',partition)
    config.spTrace_queue.set('nbundle',nbundle)
    slurmppn = int(config.spTrace_queue.get('max.ppn'))        
    config.spTrace_queue.set('ppn',min(slurmppn,max(len(mjd),4)))
    if maxjobs is None:
        maxjobs = config.spTrace_queue.get('max.jobs')
    if maxjobs is not None:
        if maxjobs < config.spTrace_queue.get('ppn'):
            config.spTrace_queue.set('ppn', maxjobs)
    config.spTrace_queue.set('no_submit', not config.pipe['Stage.run_spTrace'])


def run_spTrace():
    setup_run(mjd=config.pipe['fmjdselect.mjd'])
    obs = config.pipe['fmjdselect.obs']
    if isinstance(obs, list):
        obs = obs[0]
    queue1 = build(config.pipe['fmjdselect.mjd'], obs)
    
def build(mjd, obs):
    mjd = np.atleast_1d(mjd)
    skip_plan = config.pipe['Stage.run_spTrace_plan']
    clobber = config.pipe['Clobber.clobber_spTrace']
    debug = config.pipe['reduce.debug']
    saveraw = config.pipe['reduce.saveraw']
    daily = config.pipe['fmjdselect.trace_all_mjds']
    if obs.lower() == 'lco':
        lco = True
    else:
        lco = False
    if len(mjd) > 2:
        mjdstr = f'{np.min(mjd)}-{np.max(mjd)}'
        mjd = mjd.astype(str).tolist()
    else:
        mjd = mjd.astype(str).tolist()
        mjdstr = f'{"_".join(mjd)}'
        
    label = f"{config.pipe['general.RUN2D']}/{obs.upper()}/{mjdstr}/run_spTrace"


    queue1 = None

    if not skip_plan:
        nmjds = spplanTrace(obs = obs,mjd=mjd)
        # topdir=config.pipe['general.BOSS_SPECTRO_REDUX'],
        #                     run2d=config.pipe['general.RUN2D'],
        #                     mjd=mjd, lco=lco, mjd_plans=(not daily))
        if nmjds is None:
            print('No Valid MJDs... skipping spTrace')
            return(None, None, None)
    i = 0
    if len(mjd) > 1: nthreads = 0
    elif len(mjd) == 1: nthreads = 2
    else: nthreads = 4
    
    if not config.spTrace_queue.get('shared'):
        if int(config.spTrace_queue.get('max.ppn'))*int(config.spTrace_queue.get('nodes')) > int(config.spTrace_queue.get('ppn')):
            nthreads = np.floor((int(config.spTrace_queue.get('max.ppn'))-1)*config.spTrace_queue.get('nodes').nodes/len(mjd)).astype(int)
            if nthreads == 1:
                nthreads = 0
            
    for mj in mjd:
        if not ptt.exists(ptt.join(config.pipe['general.BOSS_SPECTRO_REDUX'],
                                   config.pipe['general.RUN2D'],'trace',f'{mj}')):
            continue
        if not daily:
            if not ptt.exists(ptt.join(config.pipe['general.BOSS_SPECTRO_REDUX'],
                                       config.pipe['general.RUN2D'],'trace',f'{mj}',f"spPlanTrace-{mj}_{obs.upper()}.par")):
                continue
        idl = f'spbuild_traceflat, plan2d="spPlanTrace-{mj}_{obs.upper()}.par"'
        if debug:
            idl = idl +', /debug'
        if saveraw:
            idl = idl +', /saveraw'

        cmd = []
        cmd.append('# Auto-generated batch file '+datetime.datetime.now().strftime("%c"))
        cmd.append("#- Echo commands to make debugging easier")
        cmd.append("set -o verbose")

        script = f"idl -e '{idl}'"
        cmd.append(script)
        cmd.append(f"boss_arcs_to_traces --mjd {mj} --obs {obs.lower()} --vers {config.pipe['general.RUN2D']} --threads {nthreads}")
        cmdfile =  ptt.join(config.pipe['general.BOSS_SPECTRO_REDUX'],
                            config.pipe['general.RUN2D'],
                            'trace',f'{mj}',f"run_spTrace_{mj}_{obs.upper()}")
        with open(cmdfile,'w') as r:
            for c in cmd:
                r.write(c+'\n')
        if i == 0:
            if config.spTrace_queue.get('nodes') > 1:
                label = label.replace('/','_')
            print(config.spTrace_queue.to_str())
            queue1 = Queue(config.spTrace_queue, verbose=True)
            queue1.create(**config.spTrace_queue.to_dict(label))
        queue1.append('source '+cmdfile,outfile = cmdfile+".o.log",
                                        errfile = cmdfile+".e.log")
        i = i+1
    if len(mjd) == 0 or queue1 is None:
        print('No Valid MJDs... skipping spTrace')
        config.spTrace_queue.set('nosubmit',True)
    if not config.spTrace_queue.get('nosubmit'):
        queue1.commit(hard=True,submit=not config.spTrace_queue.get('no_submit'))
        return(queue1, cmdfile+".o.log", cmdfile+".e.log")
    else:
        return(None, None, None)
