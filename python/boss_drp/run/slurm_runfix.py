
from boss_drp.utils.daily_log import (summary, daily_log_to_file)
from boss_drp.utils.daily_log.Flag import (stopped, Error_warn, running, NoExp)
from boss_drp.field import Field
from boss_drp import daily_dir
from boss_drp.run.uubatchpbs import make_run_cmd
from boss_drp.Config import config
from boss_drp.run.queue import Queue
import boss_drp 
from boss_drp.utils import jdate
from boss_drp.utils.splog import splog

jdate = jdate.astype(int)

from glob import glob
import numpy as np
import os.path as ptt
from os import getenv
import datetime
import itertools


def _check_step(step, cf):
    if step in ['plans','spfibermap']:
        cf = ('all',None)
    elif step == 'spreduce2D':
        cf = ('spec2d','reduce2d')
    elif step == 'specombine':
        cf = ('comb','combine')
    elif step in ['spreduce1d','spXCSAO']:
        cf = ('spec1d','analyze')
    elif step == ['Fieldlist']:
        cf = ('post','fieldlist')
    elif step == ['Fieldmerge']:
        cf = ('post','fieldmerge')
    elif step == ['Reformat']:
        cf = ('reformat','specFiles')
    elif step == 'SpCalib':
        cf = 'spcalib'
    else:
        pass
    return(cf)

_verify_step = {'reduce':'spreduce2d','coadd':'rm_combine_script','spec1d':'spreduce1d_empca'}

def build_fix(topdir,run2d,run1d, mjd, obs, directory, epoch = False,
                full=False, fix_running = False, no_write=False):
    mjds = np.atleast_1d(mjd)
    obs = np.atleast_1d(obs).tolist()
    obs = '???' if len(obs) > 1 else obs[0].upper()
    redux = []
    for mjd in mjds:
        _df = summary(directory, run2d, epoch = epoch, error=True,
                        mjd = mjd, obs =  obs, html = False)
        if _df is None:
            return(redux)
        for i, row in _df.iterrows():
            try:
                field = row['Field'].split(">")[1].split("<")[0]
            except IndexError:
                field = row['Field']
            fredux = ptt.join(Field(topdir, run2d, field, epoch=epoch).dir(),
                            f"redux-{field}-{row['MJD']}")

            tredux = ptt.join(ptt.dirname(fredux),ptt.basename(fredux.replace('redux','redux_fix')))

            if no_write:
                if ptt.exists(tredux):
                    redux.append(tredux)
                continue

            cmd = []
            if full:
                cf = ('all',None)
            else:
                cf = None
                steps = ['plans','spfibermap','spreduce2D','specombine', 'spreduce1d','spXCSAO','Fieldlist','Fieldmerge','Reformat','SpCalib']
                for step in ['plans','spfibermap','spreduce2D','specombine',
                            'spreduce1d','spXCSAO','Fieldlist','Fieldmerge',
                            'Reformat','SpCalib']:
                    flagged = [f'color:{stopped.color}', f'color:{Error_warn.color}']
                    flagged.append(f'color:{NoExp.color}')
                    if fix_running:
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
            flags = []
            flags.append(f"--clean {cf[0]}")
            if topdir != getenv('BOSS_SPECTRO_REDUX'):
                flags.append(f'--topdir {topdir}')
            if run2d != getenv('RUN2D'):
                flags.append(f'--run2d {run2d}')
            if run1d != getenv('RUN1D'):
                flags.append(f'--run2d {run1d}')
            if epoch:
                flags.append(f'--epoch')
            if full:
                flags.extend(['--reset','--remove_redux'])
            flags.append(f'--field {field}')
            flags.append(f"--mjd {row['MJD']}")
            flags = ' '.join(flags)
            clean_cmd = f"boss_drp clean run {flags}"

            with open(fredux,'r') as f:
                fullcmd = f.readlines()
                fullcmd = [x.replace('\n','') for x in fullcmd]

            template = ptt.join(ptt.dirname(boss_drp.__file__), 'etc','templates','redux.j2')
            with open(template) as template_file:
                template_lines = template_file.readlines()

            cmd = []
            for i, line in enumerate(fullcmd):
                if line.strip() == "set -o verbose":
                    cmd.append(line)
                    break
                if 'Auto-generated' in line:
                    line = '# Auto-generated batch file '+datetime.datetime.now().strftime("%c")
                cmd.append(line)
            cmd.append("")
            cmd.append("#- Clean previous failed reduction steps")
            cmd.append(clean_cmd)


            for i, line in enumerate(template_lines):
                if cf[1] is not None:
                    if line.strip() == f"{{% if {cf[1]} %}}":
                        break
                elif line.strip() == 'set -o verbose':
                    break
                elif line.strip() == 'set -e':
                    break
            if i+1 < len(template_lines):
                template_lines = template_lines[i+1:]
            template_lines = [l.strip() for l in template_lines]
            
            keep_lines = []
            for line in template_lines:
                line = line.split()
                if len(line) == 0:
                    continue
                if line[0] in ['{%-','{%']:
                    continue
                if line[0] == 'touch':
                    for s in ['spec2d','specombine','spec1d']:
                        if s in line[1]:
                            keep_lines.append(f'touch_{s}')
                            break
                    continue
                if line[0] not in ['#-', 'echo']:
                    keep_lines.append(line[0])
                    continue
                keep_lines.append(line[1])
            keep_lines.append('module')

            for line in fullcmd:
                sline = line.strip().split()
                if len(sline) == 0:
                    cmd.append(line)
                    continue
                if (line[0] == '#') and (len(sline) == 1):
                    cmd.append(line)
                    continue
                if sline[0] in ['#-']:
                    if line in cmd:
                        continue
                    cmd.append(line)
                    continue
                if sline[0] == 'touch':
                    for s in ['spec2d','specombine','spec1d']:
                        if s not in sline[1]:
                            continue
                        if f'touch_{s}' in keep_lines:
                            cmd.append(line)
                            break
                    continue
                
                if 'verify' == sline[1]:
                    if sline[2] in ['setup', 'create_log']:
                        cmd.append(line)
                        continue
                    added = False
                    vs = _verify_step.get(sline[2],sline[2])
                    vs = vs.replace("'",'').replace('"','').replace(',','')
                    for l in cmd:
                        if vs in l:
                            cmd.append(line)
                            added=True
                            break
                    if added:
                        continue
                    cmd.append(line+' --skip')
                    continue

                if sline[0] not in ['#-', 'echo']:
                     if sline[0] in keep_lines:
                         cmd.append(line)
                         continue

                if sline[1] in keep_lines:
                    cmd.append(line)   
                    continue             
    
            cmd = [key for key, _ in itertools.groupby(cmd)]

            with open(tredux,'w') as r:
                for c in cmd:
                    r.write(c+'\n')
            redux.append(tredux)
    return(redux)


def slurm_runfix(full=False,fix_running=False):
    
    topdir = config.pipe['general.BOSS_SPECTRO_REDUX']
    run2d = config.pipe['general.RUN2D']
    run1d = config.pipe['general.RUN1D']
    obs = config.pipe['fmjdselect.obs']
    custom = config.pipe['customSettings.custom_name']
    epoch = config.pipe['fmjdselect.epoch']
    mjd = config.pipe['fmjdselect.mjd']
    mjdrange = config.pipe['fmjdselect.mjdrange']

    if mjd is None:
        if custom is None:
            if mjdrange is None:
                mjd = [jdate-1, jdate]
            else:
                mjd = []
                for r in mjdrange:
                    if isinstance(r[0], (list,tuple)):
                        for rr in r:
                            mjd.extend(range(rr[0], rr[1]+1))
                    else:
                        mjd.extend(range(r[0], r[1]+1))
                mjd = np.unique(mjd).tolist()
        else:
            fd = Field(topdir, run2d, '{custom}_{obs}',
                       custom_name = custom, 
                       custom = True)
            redux = []
            for _obs in obs:
                redux.extend(glob(fd.dir().format(custom=custom, 
                                                  obs=_obs.lower())))
            mjd = [int(x.split('-')[-1]) for x in redux]


    if epoch:
        dir_ = 'epoch'
    elif config.pipe['customSettings.custom_name'] is not None:
        dir_ = custom
    else:
        dir_ = 'daily'

    
    obsstr = '_'.join(obs).lower()
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

    if config.queue.get('nodes') > 1:
        label = label.replace('/','_')

        
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
                              fix_running=fix_running))
    if len(redux_list) == 0:
        splog.info(f"No Failures for {mjd} {','.join(obs)}")
        return
    queue1 = Queue(config.queue, key=None, verbose=True)
    queue1.create(**config.queue.to_dict(label=label))
        
    for redux in redux_list:
        cmd, log, err = make_run_cmd(redux)
        cmd = f'echo {cmd}'
        queue1.append(cmd, outfile = log, errfile = err)
    queue1.commit(hard=True, submit=(not config.queue.get('no_submit')))
