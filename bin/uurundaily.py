#!/usr/bin/env python3
import argparse
from slurm import queue
from os import getenv, makedirs, popen, chdir, getcwd
import os.path as ptt
from pathlib import Path
from uubatchpbs import uubatchpbs
from dailylogger import *
from daily_log import daily_log_email
from spplan import spplan1d, spplan2d
from load_module import load_module, load_env
import slurm_readfibermap
import run_spTrace
import slurm_Summary
from pydl.pydlutils.yanny import yanny, write_table_yanny, read_table_yanny
import numpy as np
from astropy.table import Table
import logging
import datetime
import astropy.time
import time
import sys
from glob import glob
import re
from tqdm import tqdm


jdate = str(int(float(astropy.time.Time(datetime.datetime.utcnow()).jd) - 2400000.5))

try:
    daily_dir = getenv('DAILY_DIR')
except:
    daily_dir = ptt.join(getenv('HOME'),'daily')



def printAndRun(log, cmd, idlspec2d_dir):
    log.info('Running '+cmd)
    cmd.replace('spplan', idlspec2d_dir+'/pro/spec2d/spplan')
    stream = popen(cmd)
    log.info(stream.read())
    log.info('')


def read_module(mod,king=False):
    module = load_module()
    module('purge')
    module('load', mod)
    run2d = load_env('RUN2D')
    run1d = load_env('RUN1D')
    if king:
        module('switch','slurm','slurm/kingspeak-pipelines')
    boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    boss_spectro_data_N =load_env('BOSS_SPECTRO_DATA_N')
    idlspec2d_dir =load_env('IDLSPEC2D_DIR')
    return(run2d, run1d, boss_spectro_redux, boss_spectro_data_N,idlspec2d_dir)

def get_nextmjd(mod, obs, nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')):
    try:
        nextmjds = yanny(nextmjd_file)
    except:
        nextmjds = {}
        nextmjds["NEXTMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))
    obss  = np.char.upper(nextmjds["NEXTMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) == 0:
        nextmjd = jdate
    else:
        nextmjd = nextmjds["NEXTMJD"]['mjd'][indx][0]
    return(nextmjd)

def check_complete(mod, obs, flag_file = ptt.join(daily_dir,'etc','completemjd.par')):
    try:
        nextmjds = yanny(nextmjd_file)
    except:
        nextmjds = {}
        nextmjds["COMPLETEMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))
    obss  = np.char.upper(nextmjds["COMPLETEMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["COMPLETEMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) == 0:
        nextmjd = jdate
    else:
        nextmjd = nextmjds["COMPLETEMJD"]['mjd'][indx][0]
    return(nextmjd)

def increment_nextmjd(logger, mod, obs, nextmjd, nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')):
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
    
def flag_complete(logger, mod, mjd, obs, flag_file = ptt.join(daily_dir,'etc','completemjd.par')):
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
            nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par'),
            flag_file = ptt.join(daily_dir,'etc','completemjd.par'),
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
                            ptt.join(daily_dir, 'etc','emails'), None, from_domain=from_domain)
                #email(subj = 'skipping '+str(lastmjd)+' for '+mod+' obs='+obs)
            lastmjd = lastmjd - 1
    if len(mjd) == 0:
        logger.info('MJD '+str(nextmjd)+' for run2d='+run2d+' OBS='+obs+' is not here yet')
    else:
        logger.info('MJDs for run2d='+run2d+' OBS='+obs+ ' transfered: '+str(nextmjd)) 
    return mjd

def dailysummary(logger, run2d, run1d, topdir, module, mjd, epoch=False,
                 pause=300, jobname='', no_submit=False):
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


def monitor_job(logger, queue1, pause = 300, jobname = '', return_status=False):
    percomp1 = 0
    q1done=False
    while percomp1 < 100:
        if queue1 is not None and not q1done:
            if queue1.get_job_status() is None:
                status = f'Failure in slurm queue for {jobname}'
                logger.info(status)
                
            t_percomp1 = queue1.get_percent_complete() if not q1done else 100
            if t_percomp1 != percomp1:
                percomp1 = t_percomp1
                logger.info(f'{jobname} {percomp1}% complete at {datetime.datetime.today().ctime()}')
        elif not q1done:
            percomp1 = 100
            status = f'{jobname} not submitted at {datetime.datetime.today().ctime()}'
            logger.info(status)
        if percomp1 == 100 and not q1done:
            q1done=True
            status = f'Complete {jobname} '
            logger.info(status)
        else:
            time.sleep(pause)
    if return_status:
        return logger, status
    
    return logger

def build_fibermaps(logger, topdir, run2d, plan2ds, mjd, obs, clobber= False,
                   pause=300, module =None, fast=False, no_submit = False):
    setup = slurm_readfibermap.Setup()
    setup.boss_spectro_redux = topdir
    setup.run2d = run2d
    setup.alloc = load_env('SLURM_ALLOC')
    
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

    queue1 = slurm_readfibermap.build(module, plan2ds, setup, obs = obs, daily = True,
                                      mjd = mjd, clobber=clobber, no_submit = no_submit)
    if queue1 is None:
        return logger
    if not no_submit:
        pause = 60
        logger = monitor_job(logger, queue1, pause=pause, jobname='slurm_readfibermap')
    return logger
    
def build_traceflats(logger, mjd, obs, run2d, topdir, clobber=False, pause=300, fast=False,
                     skip_plan=False, no_submit = False, module = None, from_domain="chpc.utah.edu"):
                     
    mjds = ','.join(np.atleast_1d(mjd).astype(str).tolist())

    if not no_submit:
        send_email('build_traceflats '+run2d +' MJD='+mjds +' OBS='+','.join(obs),
                    ptt.join(daily_dir, 'etc','emails'), None, logger, from_domain=from_domain)

    setup = run_spTrace.Setup()
    setup.boss_spectro_redux = topdir
    setup.run2d = run2d
    setup.alloc = load_env('SLURM_ALLOC')
    setup.partition = load_env('SLURM_ALLOC')
    setup.mem_per_cpu = 7500
    setup.walltime = '20:00:00'
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
        queue1, logfile, errfile = run_spTrace.build(mjd, ob, setup, clobber=clobber, module = module,
                                   skip_plan = skip_plan, no_submit = no_submit)
        attachments.extend([logfile,errfile])
        if not no_submit:
            logger = monitor_job(logger, queue1, pause=pause, jobname='run_spTrace')
        for mj in np.atleast_1d(mjd):
            for ob in obs:
                plan = read_table_yanny(ptt.join(topdir,run2d,'trace',str(mj),f'spPlanTrace-{mj}_{ob.upper()}.par'),'SPEXP')
                plan.convert_bytestring_to_unicode()
                arcs = plan['flavor'] == 'arc'
                for ff in plan[arcs]['name'].data.flatten():
                    ff = ff.replace('.fit','.fits').replace('sdR','spTraceTab')
                    if not ptt.exists(ptt.join(topdir,run2d,'trace',str(mj),ff)):
                        status = 'Fail'
#        spTraceFlats = glob(ptt.join(topdir,run2d,'trace',str(mj),'spTraceFlat-*'))
#        spTraceArcs  = glob(ptt.join(topdir,run2d,'trace',str(mj),'spTraceArc-*'))
#        if (len(spTraceFlats) != 2) and (len(spTraceArcs) != 2):
#            status = 'Fail'
    return logger, status, attachments
    
def build_run(skip_plan, logdir, obs, mj, run2d, run1d, idlspec2d_dir, options, topdir, today,
              module, plates = False, epoch=False, build_summary = False, pause=300,
              monitor=False, noslurm=False, no_dither=False, from_domain="chpc.utah.edu",
              traceflat=False, no_prep = False, clobber = False):
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
    else:
        mjfile = ptt.join(logdir, str(mj[0])+'-'+str(mj[-1])+es+'.log')
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

    if not skip_plan:
        lco = True if obs[0].upper() == 'LCO' else False
        try:
            plans2d = spplan2d(topdir=topdir, run2d=run2d, mjd=mj, lco=lco, plates=plates,
                               splog=logger, no_dither=no_dither, returnlist=True,
                               clobber = clobber)
            spplan1d(topdir=topdir, run2d=run2d, mjd=mj, lco=lco, plates=plates,
                     daily=True, splog=logger, clobber = clobber)
        except Exception as e: # work on python 3.x
            logger.error('Failure in building spPlans: '+ str(e))
            if monitor:
                logger.removeHandler(mjconsole)
                logger.removeHandler(mjfilelog)
                mjfilelog.close()
                mjconsole.close()
                send_email('Failure '+run2d +' MJD='+str(jdate) +' OBS='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), None, logger, from_domain=from_domain)
                logger.removeHandler(rootfilelog)
                rootfilelog.close()
            exit
    else:
        logger.info('Using old spplan files')

    if not no_prep:
        logger.info('Building spFibermaps for new spplan2ds')
        logger = build_fibermaps(logger, topdir, run2d, plans2d, mj, obs, clobber= clobber, pause=pause,
                                 fast = options['fast'], module = module, no_submit = no_prep,)
    
    if traceflat:
        logger.info('Building TraceFlats for mjd')
        logger, status, spTatt = build_traceflats(logger, mj, obs, run2d, topdir,
                                                  clobber=clobber, module = module,
                                                  pause=pause, skip_plan=skip_plan,
                                                  no_submit = no_prep,
                                                  from_domain=from_domain,
                                                  fast = options['fast'])
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
                send_email('spTrace Failure '+run2d +' MJD='+str(jdate) +' OBS='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), attachments, logger, from_domain=from_domain)
                logger.removeHandler(rootfilelog)
                rootfilelog.close()
            exit()
    else:
        spTatt = [None]
    fast_msg = '_fast' if options['fast'] else ''
    
    es = '' if not epoch else ' --epoch'

    logger.info('Running uubatchpbs.py --run2d '+run2d+' --obs '+obs[0]+' --sdssv'+fast_msg+' --email'+
                     ' --topdir '+topdir+ ' --run1d '+run1d+ es + 
                     ' --mjd '+' '.join(np.asarray(mj).astype(str).tolist()))
    logger.info('')
    queue1, redux = uubatchpbs(**options, obs=obs, run2d = run2d, run1d = run1d,
                               topdir = topdir, mjd=mj, logger=logger)

    if monitor and not noslurm:
        jobname = f"{run2d} MJD={','.join(np.asarray(mj).astype(str).tolist())} OBS={','.join(obs)}"
        logger, subj = monitor_job(logger, queue1, pause=pause,
                                   jobname = 'uubatch '+jobname,
                                   return_status = True)
        if ("not submitted" in subj) or ("Failure" in subj):
            build_summary = False
        if build_summary:
            logger, attach_summ, subj_sum = dailysummary(logger, run2d, run1d, topdir, module, mj,
                                                         epoch = epoch, pause=pause,
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
            
        #send_email(subj, ptt.join(getenv('HOME'), 'daily', 'etc','emails'), attachments, logger, from_domain=from_domain)
    logger.removeHandler(rootfilelog)
    rootfilelog.close()
    return
    


def uurundaily(module, obs, mjd = None, clobber=False, fast = False, saveraw=False, skip_plan=False,
              pause=300, nosubmit=False, noslurm=False, batch=False, debug=False, nodb=False, epoch=False,
              build_summary=False, monitor=False, merge3d=False, no_dither=False, traceflat=False,
              from_domain="chpc.utah.edu", king=False, no_prep = False):
 
    run2d, run1d, topdir, boss_spectro_data, idlspec2d_dir = read_module(module, king=king)
    if not epoch:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')
    else:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd_epoch.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')

    today = datetime.datetime.today().strftime("%m%d%Y")
    logdir = ptt.join(getenv("HOME"), "daily", "logs", obs[0].upper(),run2d.upper())


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
        options = {'MWM_fluxer'     : True,
                   'map3d'          : dmap,
                   'no_reject'      : True,
                   'no_merge_spall' : True,
                   'walltime'       : '40:00:00',
                   'include_bad'    : True,
                   'xyfit'          : True,
                   'loaddesi'       : True,
                   'kingspeak'      : king,
                   'shared'         : shared,
                   'mem_per_cpu'    : mem_per_cpu,
                   'fast'           : fast,
                   'email'          : True,
                   'nosubmit'       : nosubmit,
                   'daily'          : True,
                   'module'         : module,
                   'clobber'        : clobber,
                   'saveraw'        : saveraw,
                   'debug'          : debug,
                   'no_db'          : nodb,
                   'no_write'       : noslurm,
                   'epoch'          : epoch,
                   }
        rootlogger.info('')
        if batch is True:
            mjd = np.asarray(mjd)
            plate_mjds = mjd[np.where(mjd <  59540)[0]]
            if len(plate_mjds) >0:
                build_run(skip_plan, logdir, obs, plate_mjds.tolist(), run2d, run1d, idlspec2d_dir, options,
                          topdir, today, module, plates = True, epoch=epoch, build_summary=build_summary,
                          pause=pause, monitor=monitor, noslurm=noslurm, no_dither=no_dither,
                          from_domain=from_domain, traceflat=traceflat, no_prep = no_prep, clobber = clobber)
            fps_mjds   = mjd[np.where(mjd >= 59540)[0]]
            if len(fps_mjds) > 0:
                build_run(skip_plan, logdir, obs, fps_mjds.tolist(), run2d, run1d, idlspec2d_dir, options,
                          topdir, today, module, plates = False, epoch=epoch, build_summary=build_summary,
                          pause=pause, monitor=monitor, noslurm=noslurm, no_dither=no_dither,
                          from_domain=from_domain, traceflat=traceflat, no_prep = no_prep, clobber = clobber)

        else:
            for mj in mjd:
                if mj < 59540:
                    plates = True
                else: 
                    plates = False
                build_run(skip_plan, logdir, obs, [mj], run2d, run1d, idlspec2d_dir, options,
                          topdir, today, module, pause=pause, plates = plates, epoch=epoch,
                          build_summary=build_summary, monitor=monitor, noslurm=noslurm,
                          no_dither=no_dither, from_domain=from_domain, traceflat=traceflat,
                          no_prep = no_prep, clobber = clobber)


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

if __name__ == '__main__' :
    """
    Batch process Spectro-2D and Spectro-1D reductions based upon already-built plan files
    """
    parser = argparse.ArgumentParser(
            prog=ptt.basename(sys.argv[0]),
            description='Build Exposure Log')
    
    parser.add_argument('--apo', action = 'store_true')
    parser.add_argument('--lco', action = 'store_true')
    parser.add_argument('--module', type=str, help='Module for daily run', default='bhm/master')
    parser.add_argument('--mjd', type=int, help='Manually run for a single/list of mjd (does not update nextmjd.par)', nargs='*')
    parser.add_argument('--range_mjd', type=parseNumList, help='Manually run for a range of mjds (does not update nextmjd.par)')
    parser.add_argument('--clobber', action='store_true', help="clobber uubatchpbs run")
    parser.add_argument('--fast', action='store_true', help='turn on fast user for slurm')
    parser.add_argument('--saveraw', action='store_true', help='save sdssproc outputs')
    parser.add_argument('--debug', action='store_true', help='save extraction debug files')
    parser.add_argument('--skip_plan', action='store_true', help='Skip  createing spplan files')
    parser.add_argument('--nosubmit', action='store_true', help='Skip submitting uubatch job (ideal for allowing editting of plans)')
    parser.add_argument('--noslurm', action='store_true', help='Skip creating uubatch job')
    parser.add_argument('--batch', action='store_true', help='run for multiple mjds in a single batch')
    parser.add_argument('--nodb', action='store_true', help='skip Database operations')
    parser.add_argument('--epoch', action='store_true', help='Run Epoch Coadds')
    parser.add_argument('--summary', action='store_true', help='Build Summary Files')
    parser.add_argument('--monitor', action='store_true', help='Monitors pipeline status')
    parser.add_argument('--pause', type=int, help='Pause time (s) in status updates', default=15*60)
    parser.add_argument('--merge3d', action='store_true', help='Use prototype 3D Dustmap (in merge mode)')
    parser.add_argument('--traceflat', action='store_true', help='Build and use TraceFlats')
    parser.add_argument('--no_prep', action='store_true', help='Skip building TraceFlats and spfibermaps before pipeline run')
    parser.add_argument('--king', action='store_true', help='Use Kingspeak nodes')
    parser.add_argument('--no_dither', action='store_true', help='Skip Dither Engineering Fields')

    args = parser.parse_args()
    
    if args.lco is True:
        args.obs = ['lco']
    if args.apo is True:
        args.obs = ['apo']
    if args.range_mjd is not None:
        if args.mjd is not None:
            args.mjd.extend(args.range_mjd)
        else:
            args.mjd = args.range_mjd

    uurundaily(args.module, args.obs, mjd=args.mjd, clobber=args.clobber, fast = args.fast, saveraw=args.saveraw, skip_plan=args.skip_plan, 
            nosubmit=args.nosubmit, batch=args.batch, noslurm=args.noslurm, debug=args.debug, nodb= args.nodb, epoch = args.epoch,
            build_summary = args.summary, pause = args.pause, monitor=args.monitor, merge3d=args.merge3d, no_dither=args.no_dither,
            from_domain="chpc.utah.edu", traceflat = args.traceflat, king = args.king, no_prep = args.no_prep)
