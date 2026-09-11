from boss_drp.utils.daily_log.Flag import *
from boss_drp.field import Field
from boss_drp import daily_dir
from pathlib import Path
import json
import sys
from os import getenv



def verify_step(step, **lf_args):
    mjd1d = lf_args.get('mjd1d')
    if step == 'spfibermap':
        logfile = 'spfibermap-{field}-{mjd}.log'
    if step == 'spDiag2d':
        logfile = 'spDiag2d-{field}-{mjd}.log'
    if step == 'spDiagcomb':
        if mjd1d is None:
            logfile = 'spDiagcomb-{field}-{mjd}.log'
        else:
            logfile = 'spDiagcomb-{field}-{mjd}_{mjd1d}.log'
    if step == 'spDiag1d':
        if mjd1d is None:
            logfile = '{run1d}/spDiag1d-{field}-{mjd}.log'
        else:
            logfile = '{run1d}/spDiag1d-{field}-{mjd1d}.log'
    if step == 'spXCSAO':
        if mjd1d is None:
            logfile = '{run1d}/spXCSAO-{field}-{mjd}.log'
        else:
            logfile = '{run1d}/spXCSAO-{field}-{mjd1d}.log'
    if step == 'fieldlist':
        logfile = 'fieldlist-{field}-{mjd}.log'
    if step == 'spAll':
        if mjd1d is None:
            logfile = 'spAll-{field}-{mjd}.log'
        else:
            logfile = 'spAll-{field}-{mjd1d}.log'
    if step == 'spSpec_reformat':
        if mjd1d is None:
            logfile = 'spSpec_reformat-{field}-{mjd}.log'
        else:
            logfile = 'spSpec_reformat-{field}-{mjd1d}.log'
    if step == 'spCalib_QA':
        logfile = 'spCalib_QA-{run2d}-{field}-{mjd}.log'

    logfile = logfile.format(**lf_args)



    line = -2 if step == 'spfibermap' else -1
    with open(logfile) as f:
        try:
            last_line = f.readlines()[line]
        except:
            last_line = ''
    
    if complete.line[step] in last_line:
        with open(logfile) as f:
            lines = f.readlines()
            lines.reverse()
            
            for i, line in enumerate(lines):
                for err in errors:
                    msg = err.check(i, line, step)
                    if msg is not None:
                        if step == 'spDiag2d':
                            if 'SPCALIB: b' in msg:
                                for li, ll in enumerate(lines):
                                    noerr = noerr_cal_b.check(li, ll, step)
                                    if noerr == 'No error':
                                        msg = None
                            elif 'SPCALIB: r' in msg:
                                for li, ll in enumerate(lines):
                                    noerr = noerr_cal_r.check(li, ll, step)
                                    if noerr == 'No error':
                                        msg = None
                        if msg is not None:
                            return(err.flag,msg) 
        return(NoIssues, None)

    if ((step in ['spfibermap','spXCSAO','fieldlist','spAll','run_spTrace'])
            and '.log' in logfile):
        with open(logfile) as f:
            lines = f.readlines()
            lines.reverse()
            
            for i, line in enumerate(lines):
                for err in py_err:
                    msg = err.check(i,line,step)
                    if msg is not None:
                        return(err.flag,msg)
        return(running, None)


run_log = {
    "name":"",
    "config":{
        "field":"",
        "mjd":"",
        "mjd1d":None,
        "epoch": False,
        "custom":None,
        "run2d": None,
        "run1d": None,
    },
    "steps": {},
    "status":"pending",
    "current_step":None,
    "message":None
}


def get_logjson(run2d, field, mjd, mjd1d=None, epoch=False, custom=None, var=False):

    status_dir =  Path(daily_dir) / 'Status_logs' / f'{run2d}'
    if custom:
        status_dir = status_dir / f'{custom}'
    elif epoch:
        status_dir = status_dir / 'epoch'
    else:
        status_dir = status_dir / 'daily'
    mjd_s = mjd1d if mjd1d is not None else mjd
    mjd_grp = '{:0>3d}XX'.format(int(mjd_s)//100)
    status_dir = status_dir / mjd_grp / mjd_s 
    logjson = status_dir / f'status_{field}-{mjd}.json'
    if mjd1d is not None:
        logjson = status_dir / f'status_{field}-{mjd}_{mjd1d}.json'
    return logjson

def setup_json(run2d, field, mjd, mjd1d=None, epoch=False, custom=None, run1d=None, obs=None, **kwargs):
   
    logjson = get_logjson(run2d, field, mjd, mjd1d, epoch, custom)
    logjson.parent.mkdir(parents=True, exist_ok=True)

    logtmp = logjson.with_suffix('.tmp')
    log = run_log
    complete.set(custom)

    for step_ in complete.line:
        if step_ == 'run_spTrace':
            continue
        log['steps'][step_] = {'status':'pending','message':None}

    if run1d is None:
        run1d = run2d
    name = logjson.name
    while Path(name).suffix:
        name = Path(name).stem.replace('status_','')
    log['name'] = name
    log['config']['obs'] = obs
    log['config']['field'] = field
    log['config']['mjd'] = mjd
    log['config']['mjd1d'] = mjd1d
    log['config']['epoch'] = epoch
    log['config']['custom'] = custom
    log['config']['run2d'] = run2d
    log['config']['run1d'] = run1d

    with open(logtmp, 'w', encoding="utf-8") as file:
        json.dump(log, file, indent=4)
    logtmp.rename(logjson)
    

steps_key = {"readfibermap":'spfibermap',
            "reduce":'spDiag2d',
            "coadd":'spDiagcomb',
            "spec1d":'spDiag1d',
            "run_pyxcsao":'spXCSAO',
            "fieldlist":'fieldlist',
            "fieldmerge":'spAll',
            "spspec_reformat":'spSpec_reformat',
            "spcalib_qa":"spCalib_QA",
            "all":"all"}

def verify(step, run2d, field, mjd, epoch=False, custom = None, skip=False, 
           verbose=False, mjd1d=None, **kwargs):
    step = steps_key[step.lower()]
    logjson = get_logjson(run2d, field, mjd, epoch=epoch, custom=custom, mjd1d=mjd1d)

    logtmp = logjson.with_suffix('.tmp')

    with open(logjson, "r", encoding="utf-8") as file:
        log = json.load(file)

    complete.set(log['config']['custom'])

    status_f = 0
    keys = list(log['steps'].keys())
    for step_, status in log['steps'].items():
        i = keys.index(step_)
        next_key = keys[i + 1] if i + 1 < len(keys) else None
        if (step != 'all') and (step != step_):
            if status['status'] in ['pending','failed']:
                if  status['status'] == 'failed':
                    log['status'] = 'failed'
                    status_f = 1
                break
        elif (step != step_):
            if status['status'] in ['failed']:
                log['status'] = 'failed'
                status_f = 1
                break
        if (step == step_) or (step == 'all'):
            if skip:
                if log['steps'][step_]['status'] == 'pending':
                    log['steps'][step_]['status'] = 'skipped'
                    log['current_step'] = next_key
                    if verbose: 
                        print(f'Skipping {step_} (Marking as skipped)')
                continue
            flag, msg = verify_step(step_, **log['config'])
            if flag in [stopped, NoExp]:
                log['steps'][step_]['status'] = 'failed'
                status_f = 1
            else:
                log['steps'][step_]['status'] = 'complete'
            log['steps'][step_]['message'] = msg
            if log['message'] is None:
                log['message'] = msg
            elif msg is None:
                pass
            else:
                log['message'] = ', '.join([log['message'], msg])

            log['current_step'] = next_key
            if log['steps'][step_]['status'] == 'failed':
                log['status'] = log['steps'][step_]['status']

            if status_f == 1:
                break
    if step != 'all' and verbose:
        print(f'Status of {step}: {log["steps"][step]}')
        print(f'Error code: {status_f}')
    elif step != 'all':
        if log["steps"][step]['status'] == 'failed':
            print(f'Failed Step: {step}')
    message = []
    log['status'] = 'complete'
    for step_, status in log['steps'].items():
        if status['message'] is not None:
            message.append(f'{step_}: {status["message"]}')
        if status['status'] in ['skipped','complete']:
            continue
        if status['status'] in ['pending']:
            log['status'] = 'pending'
            log['current_step'] = step_
            continue
        if status['status'] in ['failed']:
            log['status'] = 'failed'
            log['current_step'] = step_
            status_f = 1
            break
    log['message'] = ', '.join(message)

    with open(logtmp, 'w', encoding="utf-8") as file:
        json.dump(log, file, indent=4)
    logtmp.rename(logjson)
    sys.exit(status_f)

