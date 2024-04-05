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
from load_module import load_module, load_env
from pydl.pydlutils.yanny import read_table_yanny
import numpy as np
from field import field_to_string
from glob import glob



class Setup:
    def __init__(self):
        self.boss_spectro_redux = None
        self.run2d = None
        self.alloc = None
        self.nodes = 1
        self.ppn = None
        self.mem_per_cpu = None
        self.walltime = None
        self.shared = False
        self.partition
        
    def __repr__(self):
        return self.__str__()
            
    def __str__(self):
        return (f"boss_spectro_redux: {self.boss_spectro_redux} \n"    +
                f"run2d: {self.run2d} \n"    +
                f"alloc: {self.alloc} \n"    +
                f"partition: {self.partition} \n"    +
                f"nodes: {self.nodes} \n"    +
                f"ppn: {self.ppn} \n"    +
                f"mem_per_cpu: {self.mem_per_cpu} \n"    +
                f"walltime: {self.walltime} \n"+
                f"shared: {self.shared}");
        

def read_mod(mod,kingspeak=False):
    module = load_module()
    module('purge')
    module('load', mod)
    if kingspeak:
        module('switch', 'slurm', 'slurm/kingspeak-pipelines')
    setup = Setup()
    setup.boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    setup.run2d = load_env('RUN2D')
    setup.alloc = load_env('SLURM_ALLOC')
    setup.nodes = 1 #load_env('SLURM_NODES')
    setup.ppn = load_env('SLURM_PPN')
    setup.mem_per_cpu = load_env('SLURM_MEM_PER_CPU')
    setup.walltime = load_env('SLURM_WALLTIME')
    setup.shared = False if kingspeak else True
    setup.partition = load_env('SLURM_ALLOC')

    return(setup)




def slurm_readfibermap(module='bhm/master', walltime = '40:00:00', mem = 32000, kingspeak=False,
                      mjd=None, mjdstart=None, mjdend=None, obs=['apo','lco'], ppn=20,
                      clobber=False, dr19 = False):

    setup = read_mod(module, kingspeak=kingspeak)
    
    if ppn is not None:
        maxppn = min(20,int(setup.ppn))
        if ppn > int(setup.ppn):
            print(f'Resetting to {maxppn} due to cluster configuration')
            ppn = int(setup.ppn)
        if ppn > 20:
            if input('Are you sure you want to run more then 20 in parallel? (y/[n])  ').lower() == 'y':
                pass
            else:
                ppn = maxppn
                print(f'Reseting to {maxppn} in parallel')
        setup.ppn = ppn

    if walltime is not None:
        setup.walltime = walltime
    if mem is not None:
        setup.mem_per_cpu = mem
    

    daily_dir = getenv('DAILY_DIR')
    if daily_dir is None: daily_dir = ptt.join(getenv('HOME'), "daily")
    
    plan2ds = glob(ptt.join(setup.boss_spectro_redux,setup.run2d,
                            field_to_string(0).replace('0','?'),'spPlan2d*.par'))

    qu = build(module, plan2ds, setup, clobber=clobber, daily_dir=daily_dir,
                mjd=mjd, mjdstart= mjdstart, mjdend=mjdend, obs = obs, dr19 = dr19)

def build(module, plan2ds, setup, clobber=False, daily=False, dr19=False,
            mjd=None, mjdstart= None, mjdend=None, no_submit=False,
            obs = ['apo','lco'], daily_dir=ptt.join(getenv('HOME'), "daily")):
    i = 0
    title = 'readfibermap_'+setup.run2d
    #if not daily:
    log = ptt.join(daily_dir, "logs", "readfibermap", setup.run2d, "readfibermap_")
    makedirs(ptt.join(daily_dir, "logs", "readfibermap", setup.run2d), exist_ok = True)
    if not no_submit:
        print(setup)
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

        if i == 0:
            if not no_submit:
                queue1 = queue(key=None, verbose=True)
                queue1.create(label = title, nodes = setup.nodes, ppn = setup.ppn,
                     walltime = setup.walltime, alloc=setup.alloc, partition = setup.partition,
                     mem_per_cpu = setup.mem_per_cpu, shared = setup.shared)
            else:
                queue1 = None
        
#        if not daily:
        thislog = log+ptt.basename(plan2d).replace('spPlan2d-','').replace('.par','')
#        else:
#            thislog = plan2d.replace('spPlan2d-','spfibermap').replace('.par','')
        drf = '' if not dr19 else ' --dr19'
        thiscmd = (f"module purge ; module load {module} ; cd {ptt.dirname(plan2d)} ; " +
                  f"readfibermaps.py --spplan2d {ptt.basename(plan2d)} {drf}")
        if clobber:
            thiscmd = thiscmd+' --clobber'
        
        if not no_submit:
            queue1.append(thiscmd, outfile = thislog+".o.log", errfile = thislog+".e.log")
        i = i+1
    if i > 0:
        if not no_submit:
            queue1.commit(hard=True,submit=True)
    
    return(queue1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create daily field merge slurm job')
    parser.add_argument('--module', '-m', default = 'bhm/master', help = 'module file to use (ex bhm/master[default] or bhm/v6_0_9)')
    parser.add_argument('--clobber', action='store_true', help='Clobber spfibermaps')
    parser.add_argument('--apo', action='store_true', help='run apo')
    parser.add_argument('--lco', action='store_true', help='run lco')
    parser.add_argument('--dr19', action='store_true', help='Limit targeting flags to DR19 cartons')

    mjdgroup = parser.add_argument_group('Select MJDs')
    mjdgroup.add_argument('--mjd', nargs='*', help='MJD dates to reduce; default="*"', type=int)
    mjdgroup.add_argument('--mjdstart', help='Starting MJD', default=None, type=int)
    mjdgroup.add_argument('--mjdend', help='Ending MJD', default=None, type=int)

    slurmgroup = parser.add_argument_group('Slurm Options')
    slurmgroup.add_argument('--kingspeak', action='store_true', help='Use kingspeak rather then notchpeak')
    slurmgroup.add_argument('--mem_per_cpu', help='Memory allocated per CPU', type=str)
    slurmgroup.add_argument('--walltime', help='Wall time in hours', type=str)
    slurmgroup.add_argument('--ppn', help='Number of processors per node', type=int)
  
    args = parser.parse_args()
    
    obs = []
    if args.apo:
        obs.append('apo')
    if args.lco:
        obs.append('lco')
    if len(obs)== 0:
        obs = ['apo','lco']
                              
                      
    slurm_readfibermap(module=args.module, walltime = args.walltime, mem = args.mem_per_cpu,
                       kingspeak=args.kingspeak, mjd=args.mjd, mjdstart = args.mjdstart,
                       mjdend = args.mjdend, obs = obs, ppn = args.ppn, clobber =args.clobber,
                       dr19 = args.dr19)
