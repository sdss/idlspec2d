#!/usr/bin/env python3
from boss_drp.prep.spplan_trace import spplanTrace
from boss_drp.utils import load_env

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
from os import path as ptt
import numpy as np
import datetime


run2d = None
topdir = None


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
        self.partition = None
        self.bundle = None
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


def run_spTrace(mjd, obs, run2d, topdir, nodes = 1, clobber=False, alloc='sdss-np',
                partition =None, debug = False, skip_plan=False, no_submit=False,
                nbundle = None, daily=False, walltime = None):
    setup = Setup()
    setup.boss_spectro_redux = topdir
    setup.run2d = run2d
    setup.nodes = nodes
    setup.alloc = alloc
    setup.nbundle = nbundle
    if setup.nbundle is not None:
        setup.bundle = True
    if partition is None:
        setup.partition = alloc
    else:
        setup.partition = partition
    if 'sdss-kp' in setup.alloc:
        slurmppn = int(load_env('SLURM_PPN'))//2
        setup.mem_per_cpu = 3750
        setup.ppn = min(slurmppn,max(len(mjd),2))*2
    else:
        slurmppn = int(load_env('SLURM_PPN'))
        setup.mem_per_cpu = 7500
        setup.ppn = min(slurmppn,max(len(mjd),2))
    setup.shared = False if 'kp' in alloc else True
    if walltime is None:
        setup.walltime = '72:00:00'
    else:
        setup.walltime = walltime
#    setup.mem_per_cpu = 7500
    queue1 = build(mjd, obs, setup, clobber=clobber, skip_plan=skip_plan,
                   debug = debug, no_submit=no_submit, daily=daily)
    
def build(mjd, obs, setup, clobber=False, no_submit=False, skip_plan=False,
          debug = False, daily=False):
    mjd = np.atleast_1d(mjd)
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
        
    label = f'{setup.run2d}/{obs.upper()}/{mjdstr}/run_spTrace'


    queue1 = None

    if not skip_plan:
        nmjds = spplanTrace(topdir=setup.boss_spectro_redux,run2d=setup.run2d,
                            mjd=mjd, lco=lco, mjd_plans=(not daily))
        print(nmjds)
        if nmjds is None:
            print('No Valid MJDs... skipping spTrace')
            return(None, None, None)
    i = 0
    for mj in mjd:
        if not ptt.exists(ptt.join(setup.boss_spectro_redux,setup.run2d,'trace',f'{mj}')):
            continue
        idl = f'spbuild_traceflat, plan2d="spPlanTrace-{mj}_{obs.upper()}.par"'
        if debug:
            idl = idl +', /debug'
            
        cmd = []
        cmd.append('# Auto-generated batch file '+datetime.datetime.now().strftime("%c"))
        cmd.append("#- Echo commands to make debugging easier")
        cmd.append("set -o verbose")

        script = f"idl -e '{idl}'"
        cmd.append(script)
        cmd.append(f"boss_arcs_to_traces --mjd {mj} --obs {obs.lower()} --vers {setup.run2d} --threads 4")
        cmdfile =  ptt.join(setup.boss_spectro_redux,setup.run2d,'trace',f'{mj}',f"run_spTrace_{mj}_{obs.upper()}")
        with open(cmdfile,'w') as r:
            for c in cmd:
                r.write(c+'\n')
        if i == 0 and not no_submit:
            if setup.nodes > 1:
                label = label.replace('/','_')
            print(setup)
            queue1 = queue(verbose=True)
            if setup.nbundle is not None:
                setup.bundle = True
            queue1.create(label=label,nodes=setup.nodes,ppn=setup.ppn,
                          shared=setup.shared, walltime=setup.walltime,
                          alloc=setup.alloc, partition=setup.partition,
                          mem_per_cpu = setup.mem_per_cpu,
                          bundle = setup.bundle, nbundle = setup.nbundle)
        if not no_submit:
            queue1.append('source '+cmdfile,outfile = cmdfile+".o.log",
                                            errfile = cmdfile+".e.log")
        elif noslurm:
            print('source '+cmdfile)
        i = i+1
    if len(mjd) == 0 or queue1 is None:
        print('No Valid MJDs... skipping spTrace')
        no_submit = True
    if not no_submit:
        queue1.commit(hard=True,submit=True)
        return(queue1, cmdfile+".o.log", cmdfile+".e.log")
    else:
        return(None, None, None)
