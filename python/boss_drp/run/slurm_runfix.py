
from boss_drp.utils.daily_log import (summary, daily_log_to_file)
from boss_drp.utils.daily_log.Flag import (stopped, Error_warn, running, NoExp)
from boss_drp.field import field_dir
from boss_drp import daily_dir
from boss_drp.run.uubatchpbs import make_run_cmd
from boss_drp.utils import load_env

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
    warnings.warn('No slurm package installed: printing command to STDOUT/logfile for manual run',SlurmWarning)
    noslurm = True

import numpy as np
import os.path as ptt
import datetime

def _check_step(step, cf):
    if step in ['plans','spfibermap']:
        cf = 'all'
    elif step == 'spreduce2D':
        cf = 'spec2d'
    elif step == 'specombine':
        cf = 'comb'
    elif step in ['spreduce1d','spXCSAO']:
        cf = 'spec1d'
    elif step == ['Fieldlist']:
        cf = 'post'
    elif step == ['Fieldmerge']:
        cf = 'merge'
    elif step == ['Reformat']:
        cf = 'reformat'
    elif step == 'SpCalib':
        cf = 'spcalib'
    else:
        pass
    return(cf)

def build_fix(topdir,run2d,run1d, mjd, obs, directory, epoch = False,
                full=False, running = False):
    mjds = np.atleast_1d(mjd)
    obs = np.atleast_1d(obs).tolist()
    obs = '???' if len(obs) > 1 else obs[0].upper()
    redux = []
    for mjd in mjds:
        _df = summary(directory, run2d, epoch = epoch, error=True,
                        mjd = mjd, obs =  obs, html = False)
        if _df is None:
            return(redux)
        ef = ' --epoch'if epoch else ''
        for i, row in _df.iterrows():
            cmd = []
            cff = ' --reset --remove_redux' if full else ''
            if full:
                cf = 'all'
            else:
                cf = None
                for step in ['plans','spfibermap','spreduce2D','specombine',
                             'spreduce1d','spXCSAO','Fieldlist','Fieldmerge',
                             'Reformat','SpCalib']:
                    flagged = [f'color:{stopped.color}', f'color:{Error_warn.color}']
                    if running:
                        flagged.append(f'color:{running.color}')
                    
                    for f in flagged:
                        if f in row[step]:
                            cf = _check_step(step, cf)
                        if cf is not None:
                            break
                    if cf is not None:
                        break
                if cf is None:
                    continue
                            
            clean_cmd = f"clean --clean {cf} --topdir {topdir} --run2d {run2d} {ef}{cff} --field {row['Field']} --mjd {row['MJD']}"

            fredux = ptt.join(field_dir(ptt.join(topdir,run2d),row['Field']),ef,
                                                 f"redux-{row['Field']}-{row['MJD']}")
            with open(fredux,'r') as f:
                fullcmd = f.readlines()
                fullcmd = [x.replace('\n','') for x in fullcmd]
            cmd = []
            cmd.append('# Auto-generated batch file '+datetime.datetime.now().strftime("%c"))

            for line in fullcmd:
                if cmd[-1] == 'set -o verbose':
                    cmd.append('')
                    cmd.append(clean_cmd)
                    
                if ('cd ' in line or 'export' in line or
                    '#- Echo commands to make debugging easier' in line or
                    'healpix' in line or 'set -o verbose' in line):
                    cmd.append(line)
                elif len(line) == 0:
                    cmd.append(line)
                elif 'readfibermaps' in line and cf in ['all']:
                    cmd.append(line)
                elif ('spec2d' in line or 'spreduce2d' in line) and cf in ['all','spec2d']:
                    cmd.append(line)
                elif (('specombine' in line or 'spPlancombepoch' in line or 'rm_combine_script' in line) and
                      cf in ['all','spec2d','comb']):
                    cmd.append(line)
                elif ('spec1d' in line or 'spreduce1d_empca' in line or
                    'run_PyXCSAO' in line) and cf in ['all','spec2d','comb','spec1d']:
                    cmd.append(line)
                elif (('Field List' in line or 'fieldlist' in line) and
                      cf in ['all','spec2d','comb','spec1d','post']):
                    cmd.append(line)
                elif (('Field spAll' in line or 'fieldmerge' in line) and
                      cf in ['all','spec2d','comb','spec1d','post','merge']):
                    cmd.append(line)
                elif (('final spectra file' in line or 'spSpec_reformat' in line) and
                      cf in ['all','spec2d','comb','spec1d','post','merge','reformat']):
                    cmd.append(line)
                elif (('Spectro-Photometric' in line or 'spcalib_qa' in line) and
                      cf in ['all','spec2d','comb','spec1d','post','merge','reformat','spcalib']):
                    cmd.append(line)
        
            with open(ptt.join(ptt.dirname(fredux),ptt.basename(fredux.replace('redux','redux_fix'))),'w') as r:
                for c in cmd:
                    r.write(c+'\n')
            redux.append(ptt.join(ptt.dirname(fredux),ptt.basename(fredux.replace('redux','redux_fix'))))
    return(redux)


def slurm_runfix(topdir, run2d, run1d, mjd, obs,full=False,
                 epoch=False, custom=None, running=False,
                 nosubmit = False, nbundle = None, nodes = 1,
                 shared=False, walltime = '24:00:00', no_write=False):
    if noslurm:
        no_write= True
    if epoch:
        dir_ = 'epoch'
    elif custom is not None:
        dir_ = custom
    else:
        dir_ = 'daily'

    
    obsstr = '_'.join(obs).upper()
    mjdstr = str(mjd[0]) if len(mjd) == 1 else 'batch'
    if custom is None:
        if not epoch:
            label = f'{run2d}/{obsstr}/{mjdstr}/daily/fix/'
        else:
            label = f'{run2d}/{obsstr}/{mjdstr}/epoch/fix/'
    else:
        if not epoch:
            label = f'{run2d}/{obsstr}/{mjdstr}/{custom}/fix/'
        else:
            label = f'{run2d}/{obsstr}/{mjdstr}/epoch/{custom}/fix/'

    if nodes > 1:
        label = label.replace('/','_')
    
    
    ppn = load_env('SLURM_PPN')
    alloc = load_env('SLURM_ALLOC')
    partition = alloc
    mem_per_cpu = load_env('SLURM_MEM_PER_CPU')
    if shared:
        if 'sdss-kp' in alloc:
            shared = False
    if shared:
        ppn = 10
        
    directory = ptt.join(daily_dir, 'logs', 'Status', dir_, run2d)
    redux_list = []
    for mjd in mjd:
        for ob in obs:
            
            daily_log_to_file(ob, mjd, topdir=topdir, run2d=run2d,
                            run1d=run1d, redux=None, html_log=None,
                            summary=True, epoch=epoch,
                            custom = custom)
            redux_list.extend(build_fix(topdir, run2d, run1d, mjd, ob,
                              directory, epoch = epoch, full = full,
                              running=running))
    if len(redux_list) == 0:
        print(f"No Failures for {mjd} {','.join(obs)}")
        return
    
    if not no_write:
        queue1 = queue(key=None, verbose=True)
        bundle = True if nbundle is not None else False
        queue1.create(label=label, nodes=str(nodes), ppn=str(ppn),
                      partition = partition, alloc=alloc, shared=shared,
                      walltime=walltime, mem_per_cpu=mem_per_cpu,
                      nbundle = nbundle, bundle = bundle,)
    rlist = redux_list.copy()
    for i in range(len(redux_list)):
        cmd, log, err = make_run_cmd(redux_list[0])
        redux_list.pop(0)
        if not no_write:
            queue1.append(cmd, outfile = log, errfile = err)
        else:
            print(f'{cmd}  > {log} 2> {err}')
    if not no_write:
        queue1.commit(hard=True, submit=(not nosubmit))
