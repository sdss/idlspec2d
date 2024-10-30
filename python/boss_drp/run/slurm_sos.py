#!/usr/bin/env python3
from boss_drp.utils import load_env
from boss_drp.utils import jdate
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
from os import getenv, makedirs, popen, chdir, getcwd
import os.path as ptt
from datetime import date, datetime
import astropy.time
import subprocess
import io
import numpy as np
from glob import glob



class Setup:
    def __init__(self):
        self.boss_spectro_redux = None
        self.run2d = None
        self.alloc = None
        self.nodes = None
        self.ppn = None
        self.mem_per_cpu = None
        self.walltime = None
        self.shared = False
        self.flags = []
        self.bundle = False
        self.nbundle = None
        
    def __repr__(self):
        return self.__str__()
     
    def set_flags(self, no_reject = False, clobber_fibermap = False,
                no_sdssv_sn2 = False, no_arc2trace = False, forcea2t = False):
        self.flags = ['--utah', '--nodb']
        if no_reject:
            self.flags.append('--no_reject')
        if clobber_fibermap:
            self.flags.append('-f')
        if no_sdssv_sn2:
            self.flags.append('--no_sdssv_sn2')
        if no_arc2trace:
            self.flags.append('-a')
        if forcea2t:
            self.flags.append('-o')
    
    def __str__(self):
        return (f"boss_spectro_redux: {self.boss_spectro_redux} \n"    +
                f"run2d: {self.run2d} \n"    +
                f"alloc: {self.alloc} \n"    +
                f"nodes: {self.nodes} \n"    +
                f"ppn: {self.ppn} \n"    +
                f"mem_per_cpu: {self.mem_per_cpu} \n"    +
                f"walltime: {self.walltime} \n"+
                f"shared: {self.shared} \n"+
                f"bundle: {self.bundle} \n"+
                f"nbundle: {self.nbundle}");


def read_mod(nbundle = None):
    setup = Setup()
    setup.boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    setup.run2d = load_env('RUN2D')
    setup.nbundle = nbundle
    if setup.nbundle is not None:
        setup.bundle = True
    if not noslurm:
        setup.alloc = load_env('SLURM_ALLOC')
        setup.nodes = int(load_env('SLURM_NODES'))
        setup.ppn = int(load_env('SLURM_PPN'))
        setup.mem_per_cpu = load_env('SLURM_MEM_PER_CPU')
        setup.walltime = load_env('SLURM_WALLTIME')
        setup.shared = True
        setup.ntasks = int(load_env('SLURM_PPN'))//2
    else:
        setup.nodes = 1
        setup.ppn = 1
    return(setup)

def slurm_SOS(walltime = '40:00:00', mem = None, nbundle = None,
        mjd=None, mjdstart=None, mjdend=None, obs=['apo','lco'],
        ppn=20, no_submit=False, nodes=1, **flags):

    setup = read_mod(nbundle = nbundle)
    if nodes is not None:
        if nodes > setup.nodes:
            print(f'Resetting to nodes to {setup.nodes} due to cluster configuration')
        else:
            setup.nodes = nodes

    if ppn is not None:
        if ppn > int(setup.ppn):
            print(f'Resetting ppn to {setup.ppn} due to cluster configuration')
        else:
            setup.ppn = ppn

    if walltime is not None:
        setup.walltime = walltime
    if mem is not None:
        setup.mem_per_cpu = mem

    qu = build(setup, daily_dir=daily_dir, mjd=mjd, mjdstart= mjdstart,
               mjdend=mjdend, obs = obs, no_submit=no_submit, **flags)

def build(setup, mjd=None, mjdstart= None, mjdend=None, no_submit=False,
            obs = ['apo','lco'], daily_dir=daily_dir, **flags):
    if noslurm:
        no_submit = True
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
    if maxppn < setup.nodes*setup.ntasks:
        setup.ppn = maxppn
        setup.ntasks = 2 * maxppn

    if not no_submit:
        print(setup)
    setup.set_flags(**flags)

    for mj in mjd:
        for ob in obs:
            for cam in ['j']:    #'b','r']:
                if i == 0:
                    if not no_submit:
                        queue1 = queue(key=None, verbose=True)
                        if setup.nbundle is not None:
                            setup.bundle = True
                        queue1.create(label = title, nodes = setup.nodes,
                                      ppn = setup.ppn, bundle = setup.bundle,
                                      nbundle = setup.nbundle,
                                      walltime = setup.walltime, alloc=setup.alloc,
                                      cpus=setup.ntasks,shared = setup.shared,
                                      mem_per_cpu = setup.mem_per_cpu)
                    else:
                        queue1 = None

                thislog = ptt.join(log, f'SOS-{ob}-{cam}-{mj}')
                thiscmd = (f"export OBSERVATORY={ob.upper()} ; " +
                           f"SOS -{cam} -m {mj} "+" ".join(setup.flags))

                if not no_submit:
                    queue1.append(thiscmd, outfile = thislog+".o.log", errfile = thislog+".e.log")
                else:
                    print(thiscmd)
                i = i+1
    if i > 0:
        if not no_submit:
            queue1.commit(hard=True,submit=True)

    return(queue1)
