#!/usr/bin/env python3
from slurm import queue
import time
from os import path as ptt
from spplan_trace import spplanTrace
from load_module import load_module
from load_module import load_env
import argparse
import numpy as np
import spplan_trace
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
    
    def __repr__(self):
        return self.__str__()
                                                                                                        
    def __str__(self):
        return (f"boss_spectro_redux: {self.boss_spectro_redux} \n"    +
                f"run2d: {self.run2d} \n"    +
                f"alloc: {self.alloc} \n"    +
                f"nodes: {self.nodes} \n"    +
                f"ppn: {self.ppn} \n"    +
                f"mem_per_cpu: {self.mem_per_cpu} \n"    +
                f"walltime: {self.walltime} \n"+
                f"shared: {self.shared}");


def run_spTrace(mjd, obs, lco, run2d, topdir, clobber=False, alloc='sdss-np', debug = False, skip_plan=False):
    setup = Setup()
    setup.boss_spectro_redux = topdir
    setup.run2d = run2d
    setup.alloc = alloc
    if 'sdss-kp' in setup.alloc:
        slurmppn = int(load_env('SLURM_PPN'))//2
    else:
        slurmppn = int(load_env('SLURM_PPN'))
    setup.ppn = min(slurmppn,max(len(mjd),2))
    setup.shared = False if 'kp' in alloc else True
    setup.walltime = '20:00:00'
    setup.mem_per_cpu = 7500
    queue1 = build(mjd, obs, setup, clobber=False, skip_plan=skip_plan, debug = debug)
def build(mjd, obs, setup, clobber=False, no_submit=False, skip_plan=False, module=None, debug = False):
    mjd = np.atleast_1d(mjd)
    if obs.lower() == 'lco':
        lco = True
    else:
        lco = False
    if len(mjd) > 2:
        label = f'run_spTrace_{np.min(mjd)}-{np.max(mjd)}_{obs.upper()}'
        mjd = mjd.astype(str).tolist()
    else:
        mjd = mjd.astype(str).tolist()
        label = f'run_spTrace_{"_".join(mjd)}_{obs.upper()}'

    if not no_submit:
        queue1 = queue(verbose=True)
        queue1.create(label=label,nodes=setup.nodes,ppn=setup.ppn,shared=setup.shared,
                     walltime=setup.walltime,alloc=setup.alloc, mem_per_cpu = setup.mem_per_cpu)
    else:
        queue1 = None
    for mj in mjd:
        if not skip_plan:
            spplanTrace(topdir=setup.boss_spectro_redux,run2d=setup.run2d, mjd=mj, lco=lco)
        idl = f'spbuild_traceflat, plan2d="spPlanTrace-{mj}_{obs.upper()}.par"'
        if debug:
            idl = idl +', /debug'
            
        cmd = []
        cmd.append('# Auto-generated batch file '+datetime.datetime.now().strftime("%c"))
        cmd.append("#- Echo commands to make debugging easier")
        cmd.append("set -o verbose")
        if module is not None:
            cmd.append(f"module purge ; module load {module}")
            cmd.append('')
        cmd.append(f"module switch miniconda miniconda/3.9")
        cmd.append(f"module load pyvista")
        script = f"idl -e '{idl}'"
        cmd.append(script)
        cmd.append(f"boss_arcs_to_traces --mjd {mj} --obs {obs.lower()} --vers {setup.run2d}")
        print(setup)
        print(setup.boss_spectro_redux,setup.run2d,'trace',f'{mj}',f"run_spTrace_{mj}_{obs.upper()}")
        cmdfile =  ptt.join(setup.boss_spectro_redux,setup.run2d,'trace',f'{mj}',f"run_spTrace_{mj}_{obs.upper()}")
        with open(cmdfile,'w') as r:
            for c in cmd:
                r.write(c+'\n')

        if not no_submit:
            queue1.append(cmdfile,outfile = cmdfile+".o.log",
                                 errfile = cmdfile+".e.log")
    if not no_submit:
        queue1.commit(hard=True,submit=True)
    return(queue1)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mjd', type=int, required=True, nargs='*')
    parser.add_argument('--lco', action = 'store_true')
    parser.add_argument('--module', type=str, help='Module for daily run', default='work/v6_1_1-tracetweak')
    parser.add_argument('--clobber', action='store_true')
    parser.add_argument('--kings', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--skip_plan', action='store_true')
    args = parser.parse_args()

    module = load_module()
    module('purge')
    if args.kings:
        module('load', 'slurm/kingspeak-pipelines')
    else:
        module('load', 'slurm/notchpeak-pipelines')
    alloc = load_env('SLURM_ALLOC')
    module('load', args.module)
    if run2d is None:
        run2d = load_env('RUN2D')
    if topdir is None:
        topdir = load_env('BOSS_SPECTRO_REDUX')

    if args.kings:
        module('load', 'slurm/kingspeak-pipelines')
    else:
        module('load', 'slurm/notchpeak-pipelines')
    alloc = load_env('SLURM_ALLOC')
    obs = 'lco' if args.lco else 'apo'
    
    run_spTrace(args.mjd, obs, args.lco, run2d, topdir, clobber=args.clobber, alloc= alloc, debug= args.debug, skip_plan=args.skip_plan)