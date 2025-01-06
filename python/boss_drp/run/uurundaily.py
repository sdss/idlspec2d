#!/usr/bin/env python3
from boss_drp.run.uubatchpbs import uubatchpbs
from boss_drp.prep.spplan import spplan1d, spplan2d
from boss_drp.utils.daily_log import daily_log_email
from boss_drp.run import slurm_readfibermap, slurm_spTrace, slurm_Summary
from boss_drp.utils import load_env, jdate, send_email
from boss_drp.run import monitor_job
from boss_drp.field import Field
from boss_drp import daily_dir, idlspec2d_dir
from boss_drp.utils.splog import splog, Splog

import argparse
import sys
try:
    from slurm import queue
    hasslurm=True
except:
    hasslurm=False
from os import getenv, makedirs, popen,environ
import os.path as ptt
from pydl.pydlutils.yanny import yanny, write_table_yanny, read_table_yanny
import numpy as np
from astropy.table import Table
import datetime
import astropy.time
import time
from glob import glob
import re
import pandas as pd
import logging
import traceback

nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')
completemjd_file = ptt.join(daily_dir,'etc','completemjd.par')

rootlogger =  Splog(name='root')

def read_module(mem_per_cpu):
    run2d = load_env('RUN2D')
    run1d = load_env('RUN1D')
    if mem_per_cpu > 8000:
        mem_per_cpu = 8000
    shared = True
    boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    boss_spectro_data_N =load_env('BOSS_SPECTRO_DATA_N')
    return(run2d, run1d, boss_spectro_redux, boss_spectro_data_N, mem_per_cpu, shared)

def get_nextmjd(mod, obs, nextmjd_file = nextmjd_file):
    try:
        nextmjds = yanny(nextmjd_file)
    except:
        nextmjds = {}
        nextmjds["NEXTMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))
    obss  = np.char.upper(nextmjds["NEXTMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) == 0:
        nextmjd = jdate.astype(str)
    else:
        nextmjd = nextmjds["NEXTMJD"]['mjd'][indx][0]
    return(nextmjd)

def check_complete(mod, obs, flag_file = completemjd_file):
    try:
        nextmjds = yanny(nextmjd_file)
    except:
        nextmjds = {}
        nextmjds["COMPLETEMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))
    obss  = np.char.upper(nextmjds["COMPLETEMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["COMPLETEMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) == 0:
        nextmjd = jdate.astype(str)
    else:
        nextmjd = nextmjds["COMPLETEMJD"]['mjd'][indx][0]
    return(nextmjd)

def increment_nextmjd(mod, obs, nextmjd, nextmjd_file = nextmjd_file):
    try:
        nextmjds = yanny(nextmjd_file)
    except:
        nextmjds = {}
        nextmjds["NEXTMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))
    obss  = np.char.upper(nextmjds["NEXTMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) != 0: nextmjds["NEXTMJD"]['mjd'][indx] = nextmjd
    tab_nextmjds = Table(nextmjds["NEXTMJD"])
    if len(indx) == 0: tab_nextmjds.add_row([mod, nextmjd, obs])
    write_table_yanny(tab_nextmjds, nextmjd_file, tablename = "NEXTMJD", overwrite = True)
    rootlogger.info("Next MJD to wait for will be "+str(nextmjd))
    
def flag_complete(mod, mjd, obs, flag_file = completemjd_file):
    try:
        nextmjds = yanny(flag_file)
    except:
        nextmjds = {}
        pass
    if len(nextmjds) == 0:
        nextmjds = {}
        nextmjds["COMPLETEMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))
    mjd = max(np.atleast_1d(mjd).astype(int))
    obss  = np.char.upper(nextmjds["COMPLETEMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["COMPLETEMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) != 0: nextmjds["COMPLETEMJD"]['mjd'][indx] = mjd
    tab_nextmjds = Table(nextmjds["COMPLETEMJD"])
    if len(indx) == 0: tab_nextmjds.add_row([mod, mjd, obs])
    write_table_yanny(tab_nextmjds, flag_file, tablename = "COMPLETEMJD", overwrite = True)

def get_MJD(boss_spectro_data, mod, obs, run2d, epoch = False,
            nextmjd_file = nextmjd_file, flag_file = completemjd_file,
            from_domain="chpc.utah.edu"):
            
    nextmjd = get_nextmjd(mod, obs, nextmjd_file = nextmjd_file)

    if epoch:
        completemjd = check_complete(mod, obs, nextmjd_file = nextmjd_file)
        if completemjd < nextmjd:
            rootlogger.info('Daily MJD '+str(nextmjd)+' for run2d='+run2d+' OBS='+obs+' is not complete yet')
            return([])
    rootlogger.info("Looking for MJDs of or after "+str(nextmjd))
    path = ptt.join(boss_spectro_data, '?????')
    def get_key(fp):
        if not ptt.isdir(fp): return(0)
        filename = ptt.basename(fp)
        int_part = filename.split()[0]
        return(int(int_part))
    files = sorted(glob(path),key=get_key)
    lastmjd = int(ptt.basename(files[-1]))
    if epoch is True:
        mjd = [lastmjd]
    else:
        mjd = []
        while lastmjd >= int(nextmjd):
            if ptt.isdir(path.replace('?????', str(lastmjd))):
                mjd.append(lastmjd)
            else:
                send_email('skipping '+str(lastmjd)+' for '+mod+' obs='+obs,
                            ptt.join(daily_dir, 'etc','emails'), None, from_domain=from_domain)
                #email(subj = 'skipping '+str(lastmjd)+' for '+mod+' obs='+obs)
            lastmjd = lastmjd - 1
    if len(mjd) == 0:
        rootlogger.info('MJD '+str(nextmjd)+' for run2d='+run2d+' OBS='+obs+' is not here yet')
    else:
        rootlogger.info('MJDs for run2d='+run2d+' OBS='+obs+ ' transfered: '+str(nextmjd))
    return mjd

def dailysummary(run2d, run1d, topdir, module, mjd, epoch=False,
                 pause=300, jobname='', no_submit=False, obs=['apo']):
    setup = slurm_Summary.Setup()
    setup.run2d = run2d
    setup.run1d = run1d
    setup.boss_spectro_redux = topdir
    setup.epoch = epoch
    setup.custom = None
    setup.merge_only = True
    setup.backup = None
    setup.limit = None
    setup.n_iter = None
#
    queue2, title, attachments, filelog = slurm_Summary.build(setup)
    if not no_submit:
        pause = 60
        subj = monitor_job(queue1, pause=pause,
                             jobname='BOSS_Summary '+jobname,
                             return_status=True)
        rootlogger.removeHandler(filelog)
        filelog.close()
    else:
        mjd = np.atleast_1d(mjd).astype(str).tolist()
        subj = f"uubatch not submitted at {run2d} MJD={','.join(mjd)} OBS={','.join(obs)}"
    return attachments, subj


def build_fibermaps(topdir, run2d, plan2ds, mjd, obs, clobber= False,
                   pause=300, module =None, fast=False, no_submit = False,
                   nbundle = None, ):
    setup = slurm_readfibermap.Setup()
    setup.boss_spectro_redux = topdir
    setup.run2d = run2d
    setup.alloc = load_env('SLURM_ALLOC')
    setup.partition = load_env('SLURM_ALLOC')

    setup.ppn = 16
    if fast:
        setup.alloc = setup.alloc+'-fast'
    if len(plan2ds) < setup.ppn:
        setup.ppn = max([len(plan2ds), 2])
    setup.mem_per_cpu = 12000
    setup.walltime = '10:00:00'
    setup.partition = load_env('SLURM_ALLOC')
    setup.shared = True
    setup.nbundle = nbundle
    if setup.nbundle is not None:
        setup.bundle = True
    try:
        queue1 = slurm_readfibermap.build(plan2ds, setup, daily = True, obs=obs, mjd = mjd,
                                          clobber=clobber, no_submit = no_submit)
    except Exception as e:
        splog.info(traceback.format_exc())
        splog.info('Failure submitting readfibermap Jobs')
        return ('Failure submitting readfibermap Jobs')
    if queue1 is None:
        splog.info('No New Fibermaps Read')
        return (None)
    if not no_submit:
        pause = 60
        monitor_job(queue1, pause=pause, jobname='slurm_readfibermap')
    return (None)
    
def build_traceflats(mjd, obs, run2d, topdir, clobber=False, pause=300, fast=False,
                     skip_plan=False, no_submit = False, module = None, nbundle = None,
                     from_domain="chpc.utah.edu", allemail=False, **kwrds):
                     
    mjds = ','.join(np.atleast_1d(mjd).astype(str).tolist())

    if not no_submit:
        send_email('build_traceflats '+run2d +' MJD='+mjds +' OBS='+','.join(obs),
                    ptt.join(daily_dir, 'etc','emails'), None,
                    from_domain=from_domain, allemail=allemail)
    setup = slurm_spTrace.Setup()
    setup.boss_spectro_redux = topdir
    setup.run2d = run2d
    setup.alloc = load_env('SLURM_ALLOC')
    setup.partition = load_env('SLURM_ALLOC')
    setup.mem_per_cpu = 7500
    setup.walltime = '20:00:00'
    setup.nbundle = nbundle
    if setup.nbundle is not None:
        setup.bundle = True

    if fast:
        setup.alloc = setup.alloc+'-fast'
    slurmppn = int(load_env('SLURM_PPN'))
    setup.ppn = min(slurmppn,max(len(mjd),4))
    
    setup.shared = True
    
    status = 'Pass'
    attachments = []
    for ob in obs:
        queue1, logfile, errfile = slurm_spTrace.build(mjd, ob, setup, clobber=clobber,
                                                       skip_plan = skip_plan,
                                                       no_submit = no_submit,
                                                       daily=True)
        if queue1 is None:
            continue
        attachments.extend([logfile,errfile])
        if not no_submit:
            monitor_job(queue1, pause=pause, jobname='slurm_spTrace')
        for mj in np.atleast_1d(mjd):
            for ob in obs:
                plan = read_table_yanny(ptt.join(topdir,run2d,'trace',str(mj),
                                        f'spPlanTrace-{mj}_{ob.upper()}.par'),'SPEXP')
                plan.convert_bytestring_to_unicode()
                arcs = plan['flavor'] == 'arc'
                for ff in plan[arcs]['name'].data.flatten():
                    ff = ff.replace('.fit','.fits').replace('sdR','spTraceTab')
                    if not ptt.exists(ptt.join(topdir,run2d,'trace',str(mj),ff)):
                        status = 'Fail'
    return status, attachments
    
def build_run(skip_plan, logdir, obs, mj, run2d, run1d, options, topdir, today,
              plates = False, epoch=False, build_summary = False, pause=300,
              monitor=False, noslurm=False, no_dither=False, from_domain="chpc.utah.edu",
              traceflat=False, no_prep = False, clobber = False, no_fibermap = False,
              dailydir=daily_dir):
    flags = ''
    flags1d = ''
    attachments = None
    if plates is True:
        flags = flags + ', /plate_s'
        flags1d = flags1d + ', /plates'
    if obs[0].upper() == 'LCO':
        flags = flags + ', /lco'
        flags1d = flags1d + ', /lco'
    es = '' if not epoch else '_epoch'
    if len(mj) == 1:
        mjfile = ptt.join(logdir, str(mj[0])+es+'.log')
        mjsub = str(mj[0])
    else:
        mjfile = ptt.join(logdir, str(mj[0])+'-'+str(mj[-1])+es+'.log')
        mjsub = str(mj[0])+'-'+str(mj[-1])

    splog.set(name = str(mj))


    try:
        splog._log.propagate = False
    except:
        pass

    logf = 'uurundaily-'+today+'.log' if not epoch else 'uurunepoch-'+today+'.log'
    splog.add_file(ptt.join(logdir, logf))
    splog.open(logfile = mjfile)

    if skip_plan is True:
        pipeplan = False
    elif skip_plan == 'pipe':
        pipeplan = False
    else:
        pipeplan = True

    if pipeplan:
        lco = True if obs[0].upper() == 'LCO' else False
        try:
            if clobber is True:
                spPlan_clobber = True
            elif clobber == 'spPlans':
                spPlan_clobber = True
            else:
                spPlan_clobber = False
            splog.info('Creating spPlan Files')
            splog.pause_file()
            args = dict(topdir=topdir, run2d=run2d, mjd=mj, lco=lco, plates=plates,
                        no_dither=no_dither, returnlist=True,
                        clobber = spPlan_clobber, single_flat = True)
            plans2d = spplan2d(**args)
            
            args = dict(topdir=topdir, run2d=run2d, mjd=mj, lco=lco, plates=plates,
                        daily=True, clobber = spPlan_clobber, plans=plans2d)
            spplan1d(**args)
            splog.unpause_file()

        except Exception as e: # work on python 3.x
            splog.error('Failure in building spPlans: '+ str(e))
            if monitor:
                splog.close()
                send_email('Failure '+run2d +' MJD='+mjsub +' OBS='+','.join(obs),
                            ptt.join(dailydir, 'etc','emails'), mjfile, from_domain=from_domain)
                splog.close_file()
            exit
    else:
        plans2d = []
        afc = Field(topdir, run2d, '*')
        if len(mj) == 1:
            plans2d_tmp = glob(ptt.join(afc.dir(), f'spPlan2d*-{mj[0]}.par'))
        else:
            plans2d_tmp = glob(ptt.join(afc.dir(), 'spPlan2d*'))
        
        mjds = np.asarray(mj).astype(str).tolist()

        for plan2d in plans2d_tmp:
            if ptt.basename(plan2d).split('-')[-1].split('.')[0] not in mjds:
                continue
            plans2d.append(plan2d)
        splog.info('Using old spplan files')
    if len(plans2d) == 0:
        no_fibermap = True
        traceflat = False
    if not no_prep and not no_fibermap:
        if clobber is True:
            fibermap_clobber = True
        elif clobber == 'fibermap':
            fibermap_clobber = True
        else:
            fibermap_clobber = False
    
        splog.info('Building spFibermaps for spplan2ds')
        args = dict(topdir=topdir, run2d=run2d, clobber= fibermap_clobber, pause=pause,
                    fast = options['fast'], no_submit = no_prep, nbundle = options['nbundle'])
        topdir = args.pop('topdir')
        run2d  = args.pop('run2d')
        splog.pause_file()
        error = build_fibermaps(topdir, run2d, plans2d, mj, obs, **args)
        splog.unpause_file()
        if error is not None:
            splog.close()
            send_email('Failure submitting readfibermap Jobs '+mjsub+' obs='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), [mjfile], from_domain=from_domain)
            splog.close_file()
            exit()

    elif no_fibermap:
        splog.info('Skipping pre-Build of spFibermaps for spplan2ds')
    if skip_plan is True:
        Traceplan = False
    elif skip_plan == 'trace':
        Traceplan = True
    else:
        Traceplan = True
        
    if clobber is True:
        Trace_clobber = True
    elif clobber == 'trace':
        Trace_clobber = True
    else:
        Trace_clobber = False
    args = dict(active=traceflat, run2d=run2d, topdir=topdir, clobber=Trace_clobber,
                fast = options['fast'], pause=pause, skip_plan=(not Traceplan),
                no_submit = no_prep, nbundle = options['nbundle'])
    topdir = args.pop('topdir')
    run2d  = args.pop('run2d')
        
    if args['active']: #traceflat:
        splog.info('Building TraceFlats for mjd')
        splog.pause_file()
        status, spTatt = build_traceflats(mj, obs, run2d, topdir, **args)
        splog.unpause_file()
        if status == 'Fail':
            splog.error('Failure in building spTraceFlats and spTraceArcs')
            if monitor:
                attachments = []
                for f in spTatt:
                    if f is None:
                        continue
                    if ptt.exists(f):
                        attachments.append(f)
                if len(attachments) == 0:
                    attachments = None
                splog.close()
                send_email('spTrace Failure '+run2d +' MJD='+mjsub +' OBS='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), attachments, from_domain=from_domain)
                splog.close_file()
            exit()
    else:
        splog.info('Skipping Building of TraceFlats for mjd')
        spTatt = [None]
    fast_msg = '_fast' if options['fast'] else ''
    
    es = '' if not epoch else ' --epoch'

    pd_ops = pd.Series(options)
    pd_ops = pd.concat([pd.Series({'run2d':run2d,'run1d':run1d,'topdir':topdir,'epoch':epoch,
                             'mjd':np.atleast_1d(mj).tolist(), 'obs':obs}),pd_ops])

    rootlogger.console.setLevel(logging.CRITICAL + 1)
    rootlogger.info('\n'+pd_ops.to_string())
    rootlogger.console.setLevel(logging.DEBUG)
    splog.info('Running uubatchpbs --run2d '+run2d+' --obs '+obs[0]+' --sdssv'+fast_msg+' --email'+
                     ' --topdir '+topdir+ ' --run1d '+run1d+ es +
                     ' --mjd '+' '.join(np.asarray(mj).astype(str).tolist()))
    splog.info('')
    
    args = dict(**options, obs=obs, run2d = run2d, run1d = run1d, topdir = topdir,
                mjd=mj)
    
    queue1, redux = uubatchpbs(**args)
    
    if monitor and not noslurm:
        jobname = f"{run2d} MJD={','.join(np.asarray(mj).astype(str).tolist())} OBS={','.join(obs)}"
        if not options['nosubmit']:
            subj = monitor_job(queue1, pause=pause,
                                       jobname = 'uubatch '+jobname,
                                       return_status = True)
        else:
            subj = f'{jobname} not submitted at {datetime.datetime.today().ctime()}'
        if ("not submitted" in subj) or ("Failure" in subj):
            build_summary = False
        if build_summary:
            attach_summ, subj_sum = dailysummary(run2d, run1d, topdir, module, mj,
                                                         epoch = epoch, pause=pause, obs=obs,
                                                         jobname = jobname, no_submit=options['nosubmit'])
            if attach_summ is None:
                attach_summ = [None]
        else:
            attach_summ = [None]
            subj_sum = None
    else:
        attach_summ = [None]
        subj_sum = None
    splog.close()
    
    if monitor and not noslurm:
        if attachments is not None: attachments.append(mjfile)
        else: attachments = [mjfile]
        for f in spTatt:
            if f is None:
                continue
            if ptt.exists(f):
                attachments.append(f)
        for f in attach_summ:
            if f is None:
                continue
            if ptt.exists(f):
                attachments.append(f)
        if subj_sum is not None:
            subj = subj + '; ' + subj_sum
        
        for mjd in mj:
            daily_log_email(subj, attachments, obs, mjd,
                            email_file = ptt.join(daily_dir, 'etc','emails'),
                            topdir=topdir, run2d=run2d, run1d=run1d,
                            from_domain=from_domain,  redux = redux)
            
    splog.close_file()
    return
    


def uurundaily(module, obs, mjd = None, clobber=False, fast = False, saveraw=False,
              skip_plan=False, pause=300, nosubmit=False, noslurm=False, batch=False,
              debug=False, nodb=False, epoch=False, build_summary=False, monitor=False,
              merge3d=False, no_dither=False, traceflat=False, email=True,
              from_domain="chpc.utah.edu", no_prep = False, walltime='40:00:00',
              mem_per_cpu=8000, allemail=False, nbundle = None, nodist=False,
              no_fibermap=False, no_healpix=False):
 
    if not hasslurm:
        raise(Exception('No slurm package'))
    run2d, run1d, topdir, boss_spectro_data, mem_per_cpu, shared = read_module(mem_per_cpu)
    

    if not epoch:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')
    else:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd_epoch.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')

    today = datetime.datetime.today().strftime("%m%d%Y")
    logdir = ptt.join(daily_dir, "logs", obs[0].upper(),run2d.upper())


    makedirs(logdir,exist_ok=True)

    logf = 'uurundaily-'+today+'.log' if not epoch else 'uurunepoch-'+today+'.log'
    rootlogger.open(logfile = ptt.join(logdir, logf), append=True)

    
    rootlogger.debug('========================================')
    rootlogger.debug('Starting at '+datetime.datetime.today().ctime())


    if obs[0].lower() != 'apo': 
        boss_spectro_data = boss_spectro_data.replace('apo', obs[0].lower())
   
    if mjd is not None:
        manual=True
    else:
        mjd = get_MJD(boss_spectro_data, module, obs[0].upper(), run2d,
                      nextmjd_file = nextmjd_file, flag_file = flag_file,
                      epoch=epoch, from_domain=from_domain)
        manual=False
    if len(mjd) > 0:
        if manual is False:
            increment_nextmjd(module, obs[0].upper(), max(mjd)+1,
                              nextmjd_file = nextmjd_file)
        dmap = 'bayestar15' if not merge3d else 'merge3d'
        shared = True
        mem_per_cpu = 8000
        pipe_clobber = True if clobber is not False else False
        options = {'MWM_fluxer'     : True,
                   'map3d'          : dmap,
                   'no_reject'      : True,
                   'no_merge_spall' : True,
                   'walltime'       : '40:00:00',
                   'shared'         : shared,
                   'mem_per_cpu'    : mem_per_cpu,
                   'fast'           : fast,
                   'email'          : True,
                   'allemail'       : allemail,
                   'nosubmit'       : nosubmit,
                   'daily'          : True,
                   'clobber'        : pipe_clobber,
                   'saveraw'        : saveraw,
                   'debug'          : debug,
                   'no_db'          : nodb,
                   'no_write'       : noslurm,
                   'epoch'          : epoch,
                   'nbundle'        : nbundle,
                   'no_healpix'     : no_healpix,
                   'nodist'         : nodist,
                   }
        rootlogger.info('')
        if batch is True:
            mjd = np.asarray(mjd)
            plate_mjds = mjd[np.where(mjd <  59540)[0]]
            if len(plate_mjds) >0:
                build_run(skip_plan, logdir, obs, plate_mjds.tolist(), run2d, run1d,
                          options, topdir, today, plates = True,
                          epoch=epoch, build_summary=build_summary,
                          pause=pause, monitor=monitor, noslurm=noslurm, no_dither=no_dither,
                          from_domain=from_domain, traceflat=False, no_prep = no_prep,
                          clobber = clobber, dailydir = daily_dir, no_fibermap = no_fibermap)
            fps_mjds   = mjd[np.where(mjd >= 59540)[0]]
            if len(fps_mjds) > 0:
                build_run(skip_plan, logdir, obs, fps_mjds.tolist(), run2d, run1d,
                          options, topdir, today, plates = False,
                          epoch=epoch, build_summary=build_summary, pause=pause,
                          monitor=monitor, noslurm=noslurm, no_dither=no_dither,
                          from_domain=from_domain, traceflat=traceflat,
                          no_prep = no_prep, clobber = clobber,
                          dailydir = daily_dir,no_fibermap = no_fibermap)

        else:
            for mj in mjd:
                if mj < 59540:
                    plates = True
                    rtf= False
                else:
                    plates = False
                    rtf = traceflat
                build_run(skip_plan, logdir, obs, [mj], run2d, run1d, options,
                          topdir, today, pause=pause, plates = plates, epoch=epoch,
                          build_summary=build_summary, monitor=monitor, noslurm=noslurm,
                          no_dither=no_dither, from_domain=from_domain, traceflat=rtf,
                          no_prep = no_prep, clobber = clobber, dailydir = daily_dir,
                          no_fibermap = no_fibermap)

        if (not manual) and monitor:
            flag_complete(module, mjd, obs[0].upper(), flag_file = flag_file)
    rootlogger.debug('Completed at '+datetime.datetime.today().ctime())


def parseNumList(string):
    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start), int(end)+1))
