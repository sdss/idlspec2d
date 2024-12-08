#!/usr/bin/env python3
from boss_drp.utils import load_env
from boss_drp.field import field_to_string
from boss_drp.field import field_dir as get_field_dir
from boss_drp import daily_dir

import sys
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
import argparse
from os import getenv, makedirs
import os.path as ptt
from pydl.pydlutils.yanny import read_table_yanny
from glob import glob


class Setup:
    def __init__(self):
        self.boss_spectro_redux = None
        self.run2d = None
        self.alloc = None
        self.partition = None
        self.nodes = 1
        self.ppn = None
        self.mem_per_cpu = None
        self.walltime = None
        self.shared = False
        self.bundle = False
        self.nbundle = None
        
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
                f"shared: {self.shared} \n"+
                f"bundle: {self.bundle} \n"+
                f"nbundle: {self.nbundle}");
        

def setup_run(run2d=None, boss_spectro_redux=None, nbundle = None):
    setup = Setup()
    if boss_spectro_redux is None:
        setup.boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    else:
        setup.boss_spectro_redux = boss_spectro_redux
    if run2d is None:
        setup.run2d = load_env('RUN2D')
    else:
        setup.run2d = run2d
    setup.nbundle = nbundle
    if nbundle is not None:
        setup.bundle = True
    if not noslurm:
        setup.alloc = load_env('SLURM_ALLOC')
        setup.partition = load_env('SLURM_ALLOC')
        setup.nodes = 1 #load_env('SLURM_NODES')
        setup.ppn = load_env('SLURM_PPN')
        setup.mem_per_cpu = load_env('SLURM_MEM_PER_CPU')
        setup.walltime = load_env('SLURM_WALLTIME')
        setup.shared = True
    else:
        setup.ppn = 1
    return(setup)




def slurm_readfibermap(run2d=None, boss_spectro_redux=None, walltime = '40:00:00', mem = 32000,
                      mjd=None, mjdstart=None, mjdend=None, obs=['apo','lco'], ppn=20,
                      clobber=False, dr19 = False, nbundle= None):

    setup = setup_run(run2d=run2d, boss_spectro_redux=boss_spectro_redux,
                      nbundle=nbundle)
    
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
    
    
    fdir = get_field_dir(ptt.join(setup.boss_spectro_redux,setup.run2d), '*')
    plan2ds = glob(ptt.join(fdir,'spPlan2d*.par'))
    
    qu = build(plan2ds, setup, clobber=clobber, daily_dir=daily_dir,
                mjd=mjd, mjdstart= mjdstart, mjdend=mjdend, obs = obs, dr19 = dr19)

def build(plan2ds, setup, clobber=False, daily=False, dr19=False,
            mjd=None, mjdstart= None, mjdend=None, no_submit=False,
            obs = ['apo','lco'], daily_dir=daily_dir):
    i = 0
    title = 'readfibermap_'+setup.run2d
    #if not daily:
    log = ptt.join(daily_dir, "logs", "readfibermap", setup.run2d, "readfibermap_")
    makedirs(ptt.join(daily_dir, "logs", "readfibermap", setup.run2d), exist_ok = True)
    if noslurm:
        no_submit = True
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
                print(setup)
                queue1 = queue(key=None, verbose=True)
                if setup.nbundle is not None:
                    setup.bundle = True
                queue1.create(label = title, nodes = setup.nodes, ppn = setup.ppn,
                     walltime = setup.walltime, alloc=setup.alloc,
                     partition=setup.partition,
                     bundle = setup.bundle, nbundle=setup.nbundle,
                     mem_per_cpu = setup.mem_per_cpu, shared = setup.shared)
            else:
                queue1 = None
        
        thislog = log+ptt.basename(plan2d).replace('spPlan2d-','').replace('.par','')
        drf = '' if not dr19 else ' --dr19'
        thiscmd = (f"cd {ptt.dirname(plan2d)} ; " +
                  f"readfibermaps --spplan2d {ptt.basename(plan2d)} {drf}")
                
        if clobber:
            thiscmd = thiscmd+' --clobber'
        
        if not no_submit:
            queue1.append(thiscmd, outfile = thislog+".o.log", errfile = thislog+".e.log")
        elif noslurm:
            print(thiscmd)
        i = i+1
    if i > 0:
        if not no_submit:
            queue1.commit(hard=True,submit=True)
    
        return(queue1)
    else:
        print('No fibermaps Built')
    return(None)

