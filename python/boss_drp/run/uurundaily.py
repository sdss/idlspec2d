#!/usr/bin/env python3
from boss_drp.run.uubatchpbs import uubatchpbs
from boss_drp.prep.spplan import spplan1d, spplan2d
from boss_drp.utils.dailylogger import Formatter, send_email
from boss_drp.utils.daily_log import daily_log_email
from boss_drp.run import slurm_readfibermap, slurm_spTrace, slurm_Summary
from boss_drp.utils import load_env, jdate
from boss_drp.run import monitor_job
from boss_drp.field import field_dir
from boss_drp import daily_dir, idlspec2d_dir

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
import logging
import datetime
import astropy.time
import time
from glob import glob
import re
import traceback

nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')
completemjd_file = ptt.join(daily_dir,'etc','completemjd.par')


def read_module(mem_per_cpu):
    run2d = load_env('RUN2D')
    run1d = load_env('RUN1D')
    if load_env('SLURM_ALLOC').lower() == 'sdss-kp':
        king = True
        if mem_per_cpu > 3750:
            mem_per_cpu = 3750
        shared = False
    else:
        king = False
        if mem_per_cpu > 8000:
            mem_per_cpu = 8000
        shared = True
    boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    boss_spectro_data_N =load_env('BOSS_SPECTRO_DATA_N')
    return(run2d, run1d, boss_spectro_redux, boss_spectro_data_N, king, mem_per_cpu, shared)

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

def increment_nextmjd(logger, mod, obs, nextmjd, nextmjd_file = nextmjd_file):
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
    logger.info("Next MJD to wait for will be "+str(nextmjd))
    
def flag_complete(logger, mod, mjd, obs, flag_file = completemjd_file):
    try:
        nextmjds = yanny(flag_file)
    except:
        nextmjds = {}
        pass
    if len(nextmjds) == 0:
        nextmjds = {}
        nextmjds["COMPLETEMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))

    obss  = np.char.upper(nextmjds["COMPLETEMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["COMPLETEMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) != 0: nextmjds["COMPLETEMJD"]['mjd'][indx] = mjd
    tab_nextmjds = Table(nextmjds["COMPLETEMJD"])
    if len(indx) == 0: tab_nextmjds.add_row([mod, mjd, obs])
    write_table_yanny(tab_nextmjds, flag_file, tablename = "COMPLETEMJD", overwrite = True)

def get_MJD(logger, boss_spectro_data, mod, obs, run2d, epoch = False,
            nextmjd_file = nextmjd_file, flag_file = completemjd_file,
            from_domain="chpc.utah.edu"):
            
    nextmjd = get_nextmjd(mod, obs, nextmjd_file = nextmjd_file)

    if epoch:
        completemjd = check_complete(mod, obs, nextmjd_file = nextmjd_file)
        if completemjd < nextmjd:
            logger.info('Daily MJD '+str(nextmjd)+' for run2d='+run2d+' OBS='+obs+' is not complete yet')
            return([])
    logger.info("Looking for MJDs of or after "+str(nextmjd))
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
                            ptt.join(daily_dir, 'etc','emails'), None, logger, from_domain=from_domain)
                #email(subj = 'skipping '+str(lastmjd)+' for '+mod+' obs='+obs)
            lastmjd = lastmjd - 1
    if len(mjd) == 0:
        logger.info('MJD '+str(nextmjd)+' for run2d='+run2d+' OBS='+obs+' is not here yet')
    else:
        logger.info('MJDs for run2d='+run2d+' OBS='+obs+ ' transfered: '+str(nextmjd)) 
    return mjd

def dailysummary(logger, run2d, run1d, topdir, module, mjd, epoch=False,
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
    queue2, title, attachments, logger, filelog = slurm_Summary.build(setup, logger)
    if not no_submit:
        pause = 60
        logger, subj = monitor_job(logger, queue1, pause=pause,
                             jobname='BOSS_Summary '+jobname,
                             return_status=True)
        logger.removeHandler(filelog)
        filelog.close()
    else:
        mjd = np.atleast_1d(mjd).astype(str).tolist()
        subj = f"uubatch not submitted at {run2d} MJD={','.join(mjd)} OBS={','.join(obs)}"
    return logger, attachments, subj


def build_fibermaps(logger, topdir, run2d, plan2ds, mjd, obs, clobber= False,
                   pause=300, module =None, fast=False, no_submit = False,
                   nbundle = None, ):
    setup = slurm_readfibermap.Setup()
    setup.boss_spectro_redux = topdir
    setup.run2d = run2d
    setup.alloc = load_env('SLURM_ALLOC')
    setup.partition = load_env('SLURM_ALLOC')

    if 'sdss-kp' in setup.alloc:
        setup.ppn = 4
    else:
        setup.ppn = 16
        if fast:
            setup.alloc = setup.alloc+'-fast'
    if len(plan2ds) < setup.ppn:
        setup.ppn = max([len(plan2ds), 2])
    setup.mem_per_cpu = 12000
    setup.walltime = '10:00:00'
    setup.partition = load_env('SLURM_ALLOC')
    setup.shared = False if 'sdss-kp' in setup.alloc else True
    setup.nbundle = nbundle
    if setup.nbundle is not None:
        setup.bundle = True
    try:
        queue1 = slurm_readfibermap.build(plan2ds, setup, daily = True, obs=obs, mjd = mjd,
                                          clobber=clobber, no_submit = no_submit)
    except Exception as e:
        logger.info(traceback.format_exc())
        logger.info('Failure submitting readfibermap Jobs')
        return (logger, 'Failure submitting readfibermap Jobs')
    if queue1 is None:
        logger.info('No New Fibermaps Read')
        return (logger, None)
    if not no_submit:
        pause = 60
        logger = monitor_job(logger, queue1, pause=pause, jobname='slurm_readfibermap')
    return (logger, None)
    
def build_traceflats(logger, mjd, obs, run2d, topdir, clobber=False, pause=300, fast=False,
                     skip_plan=False, no_submit = False, module = None, nbundle = None,
                     from_domain="chpc.utah.edu", allemail=False, **kwrds):
                     
    mjds = ','.join(np.atleast_1d(mjd).astype(str).tolist())

    if not no_submit:
        send_email('build_traceflats '+run2d +' MJD='+mjds +' OBS='+','.join(obs),
                    ptt.join(daily_dir, 'etc','emails'), None, logger,
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
    if 'sdss-kp' in setup.alloc:
        slurmppn = int(load_env('SLURM_PPN'))//2
    else:
        if fast:
            setup.alloc = setup.alloc+'-fast'
        slurmppn = int(load_env('SLURM_PPN'))
    setup.ppn = min(slurmppn,max(len(mjd),2))
    
    setup.shared = False if 'sdss-kp' in setup.alloc else True
    
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
            logger = monitor_job(logger, queue1, pause=pause, jobname='slurm_spTrace')
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
    return logger, status, attachments
    
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
    mjfilelog = logging.FileHandler(mjfile)
    mjfilelog.setLevel(logging.DEBUG)
    mjfilelog.setFormatter(Formatter())
    mjconsole = logging.StreamHandler()
    mjconsole.setLevel(logging.DEBUG)
    mjconsole.setFormatter(Formatter())

    logf = 'uurundaily-'+today+'.log' if not epoch else 'uurunepoch-'+today+'.log'
    rootfilelog = logging.FileHandler(ptt.join(logdir, logf))
    rootfilelog.setLevel(logging.DEBUG)
    rootfilelog.setFormatter(Formatter())
    logger = logging.getLogger(str(mj))
    try:
        logger.propagate = False
    except:
        pass
    logger.addHandler(rootfilelog)
    logger.addHandler(mjfilelog)
    logger.addHandler(mjconsole)
    logger.setLevel(logging.DEBUG)
    
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
        
            args = dict(topdir=topdir, run2d=run2d, mjd=mj, lco=lco, plates=plates,
                        splog=logger, no_dither=no_dither, returnlist=True,
                        clobber = spPlan_clobber)
            plans2d = spplan2d(**args)
            
            args = dict(topdir=topdir, run2d=run2d, mjd=mj, lco=lco, plates=plates,
                        daily=True, splog=logger, clobber = spPlan_clobber, plans=plans2d)
            spplan1d(**args)
        except Exception as e: # work on python 3.x
            logger.error('Failure in building spPlans: '+ str(e))
            if monitor:
                logger.removeHandler(mjconsole)
                logger.removeHandler(mjfilelog)
                mjfilelog.close()
                mjconsole.close()
                send_email('Failure '+run2d +' MJD='+mjsub +' OBS='+','.join(obs),
                            ptt.join(dailydir, 'etc','emails'), mjfile, logger, from_domain=from_domain)
                logger.removeHandler(rootfilelog)
                rootfilelog.close()
            exit
    else:
        plans2d = []
        if len(mj) == 1:
            plans2d_tmp = glob(ptt.join(field_dir(ptt.join(topdir,run2d),'*'), f'spPlan2d*-{mj[0]}.par'))
        else:
            plans2d_tmp = glob(ptt.join(field_dir(ptt.join(topdir,run2d),'*'), 'spPlan2d*'))
        
        mjds = np.asarray(mj).astype(str).tolist()

        for plan2d in plans2d_tmp:
            if ptt.basename(plan2d).split('-')[-1].split('.')[0] not in mjds:
                continue
            plans2d.append(plan2d)
        logger.info('Using old spplan files')
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
    
        logger.info('Building spFibermaps for spplan2ds')
        args = dict(topdir=topdir, run2d=run2d, clobber= fibermap_clobber, pause=pause,
                    fast = options['fast'], no_submit = no_prep, nbundle = options['nbundle'])
        topdir = args.pop('topdir')
        run2d  = args.pop('run2d')
        logger, error = build_fibermaps(logger, topdir, run2d, plans2d, mj, obs, **args)
        if error is not None:
            logger.removeHandler(mjconsole)
            logger.removeHandler(mjfilelog)
            mjfilelog.close()
            mjconsole.close()
            send_email('Failure submitting readfibermap Jobs '+mjsub+' obs='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), [mjfile], logger, from_domain=from_domain)
            logger.removeHandler(rootfilelog)
                rootfilelog.close()
            exit()

    elif no_fibermap:
        logger.info('Skipping pre-Build of spFibermaps for spplan2ds')
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
        logger.info('Building TraceFlats for mjd')

        logger, status, spTatt = build_traceflats(logger, mj, obs, run2d, topdir, **args)
        if status == 'Fail':
            logger.error('Failure in building spTraceFlats and spTraceArcs')
            if monitor:
                attachments = []
                for f in spTatt:
                    if f is None:
                        continue
                    if ptt.exists(f):
                        attachments.append(f)
                if len(attachments) == 0:
                    attachments = None
                logger.removeHandler(mjconsole)
                logger.removeHandler(mjfilelog)
                mjfilelog.close()
                mjconsole.close()
                send_email('spTrace Failure '+run2d +' MJD='+mjsub +' OBS='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), attachments, logger, from_domain=from_domain)
                logger.removeHandler(rootfilelog)
                rootfilelog.close()
            exit()
    else:
        logger.info('Skipping Building of TraceFlats for mjd')
        spTatt = [None]
    fast_msg = '_fast' if options['fast'] else ''
    
    es = '' if not epoch else ' --epoch'

    logger.info('Running uubatchpbs --run2d '+run2d+' --obs '+obs[0]+' --sdssv'+fast_msg+' --email'+
                     ' --topdir '+topdir+ ' --run1d '+run1d+ es +
                     ' --mjd '+' '.join(np.asarray(mj).astype(str).tolist()))
    logger.info('')
    
    args = dict(**options, obs=obs, run2d = run2d, run1d = run1d, topdir = topdir,
                mjd=mj, logger=logger)
    
    queue1, redux = uubatchpbs(**args)
    
    if monitor and not noslurm:
        jobname = f"{run2d} MJD={','.join(np.asarray(mj).astype(str).tolist())} OBS={','.join(obs)}"
        logger, subj = monitor_job(logger, queue1, pause=pause,
                                   jobname = 'uubatch '+jobname,
                                   return_status = True)
        if ("not submitted" in subj) or ("Failure" in subj):
            build_summary = False
        if build_summary:
            logger, attach_summ, subj_sum = dailysummary(logger, run2d, run1d, topdir, module, mj,
                                                         epoch = epoch, pause=pause, obs=obs,
                                                         jobname = jobname, no_submit=no_submit)
            if attach_summ is None:
                attach_summ = [None]
        else:
            attach_summ = [None]
            subj_sum = None
    else:
        attach_summ = [None]
        subj_sum = None
    logger.removeHandler(mjconsole)
    logger.removeHandler(mjfilelog)
    mjfilelog.close()
    mjconsole.close()
    
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
            daily_log_email(subj, attachments, logger, obs, mjd,
                            email_file = ptt.join(daily_dir, 'etc','emails'),
                            topdir=topdir, run2d=run2d, run1d=run1d,
                            from_domain=from_domain,  redux = redux)
            
    logger.removeHandler(rootfilelog)
    rootfilelog.close()
    return
    


def uurundaily(module, obs, mjd = None, clobber=False, fast = False, saveraw=False,
              skip_plan=False, pause=300, nosubmit=False, noslurm=False, batch=False,
              debug=False, nodb=False, epoch=False, build_summary=False, monitor=False,
              merge3d=False, no_dither=False, traceflat=False, email=True,
              from_domain="chpc.utah.edu", no_prep = False, walltime='40:00:00',
              mem_per_cpu=8000, allemail=False, nbundle = None,
              no_fibermap=False, no_healpix=False):
 
    if not hasslurm:
        raise(Exception('No slurm package'))
    run2d, run1d, topdir, boss_spectro_data, king, mem_per_cpu, shared = read_module(mem_per_cpu)
    

    if not epoch:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')
    else:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd_epoch.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')

    today = datetime.datetime.today().strftime("%m%d%Y")
    logdir = ptt.join(daily_dir, "logs", obs[0].upper(),run2d.upper())


    makedirs(logdir,exist_ok=True)
    rootlogger = logging.getLogger('root')

    logf = 'uurundaily-'+today+'.log' if not epoch else 'uurunepoch-'+today+'.log'
    filelog = logging.FileHandler(ptt.join(logdir, logf))
    filelog.setLevel(logging.DEBUG)
    filelog.setFormatter(Formatter())
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(Formatter())

    rootlogger.addHandler(filelog)
    rootlogger.addHandler(console)
    rootlogger.setLevel(logging.DEBUG)

    
    rootlogger.debug('========================================')
    rootlogger.debug('Starting at '+datetime.datetime.today().ctime())


    if obs[0].lower() != 'apo': 
        boss_spectro_data = boss_spectro_data.replace('apo', obs[0].lower())
   
    if mjd is not None:
        manual=True
    else:
        mjd = get_MJD(rootlogger, boss_spectro_data, module, obs[0].upper(), run2d,
                      nextmjd_file = nextmjd_file, flag_file = flag_file,
                      epoch=epoch, from_domain=from_domain)
        manual=False
    if len(mjd) > 0:
        if manual is False:
            increment_nextmjd(rootlogger, module, obs[0].upper(), max(mjd)+1,
                              nextmjd_file = nextmjd_file)
        dmap = 'bayestar15' if not merge3d else 'merge3d'
        shared = True if not king else False
        mem_per_cpu = 8000 if not king else 3750
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
            flag_complete(rootlogger, module, mjd, obs[0].upper(), flag_file = flag_file)
    rootlogger.debug('Completed at '+datetime.datetime.today().ctime())


def parseNumList(string):
    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start), int(end)+1))
