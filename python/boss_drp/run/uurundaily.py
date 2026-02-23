#!/usr/bin/env python3
from boss_drp.run.uubatchpbs import uubatchpbs
from boss_drp.prep.spplan import spplan1d, spplan2d
from boss_drp.utils.daily_log import daily_log_email
from boss_drp.run import slurm_readfibermap, slurm_spTrace, slurm_Summary
from boss_drp.utils import load_env, jdate, send_email
from boss_drp.field import Field
from boss_drp import daily_dir, idlspec2d_dir
from boss_drp.utils.splog import splog, Splog
from boss_drp.Config import config
from boss_drp.run.queue import Queue, hasslurm
import argparse
import sys

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
import argparse

nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')
completemjd_file = ptt.join(daily_dir,'etc','completemjd.par')

rootlogger =  Splog(name='root')

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
            nextmjd_file = nextmjd_file, flag_file = completemjd_file):
            
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
                            ptt.join(daily_dir, 'etc','emails'), None)
                #email(subj = 'skipping '+str(lastmjd)+' for '+mod+' obs='+obs)
            lastmjd = lastmjd - 1
    if len(mjd) == 0:
        rootlogger.info('MJD '+str(nextmjd)+' for run2d='+run2d+' OBS='+obs+' is not here yet')
    else:
        rootlogger.info('MJDs for run2d='+run2d+' OBS='+obs+ ' transfered: '+str(nextmjd))
    return mjd

def dailysummary(module, mjd, epoch=False,
                 pause=300, jobname='', no_submit=False, obs=['apo']):
    slurm_Summary.Setup()
    config.pipe['fmjdselect.epoch'] = epoch
    config.pipe['fmjdselect.custom'] = None
    config.Summary_queue.set('no_submit', no_submit)
    config.pipe['Summary.merge_only'] = True
    config.pipe['Summary.batchwise.backup'] = None
    config.pipe['Summary.batchwise.limit'] = None
    config.pipe['Summary.batchwise.n_iter'] = None
#
    queue2, _, attachments, filelog = slurm_Summary.build()
    
    if not config.Summary_queue.get('no_submit'):
        subj = queue2.monitor_job(pause=60,jobname='BOSS_Summary '+jobname, return_status=True)

        rootlogger.removeHandler(filelog)
        filelog.close()
    else:
        mjd = np.atleast_1d(mjd).astype(str).tolist()
        subj = f"uubatch not submitted at {config.pipe['general.RUN2D']} MJD={','.join(mjd)} OBS={','.join(obs)}"
    return attachments, subj


def build_fibermaps( plan2ds, mjd, obs):

    slurm_readfibermap.setup_run()
    config.readfibermap_queue.set('ppn', 16)
    if len(plan2ds) < config.readfibermap_queue.get('ppn'):
        config.readfibermap_queue.set('ppn', max([len(plan2ds), 2]))
    config.readfibermap_queue.set('mem_per_cpu', 12000)
    config.readfibermap_queue.set('wall', '10:00:00')
    config.readfibermap_queue.set('no_submit',False)
    try:
        queue1 = slurm_readfibermap.build(plan2ds, daily = True, obs=obs, mjd = mjd)
    except Exception as e:
        splog.info(traceback.format_exc())
        splog.info('Failure submitting readfibermap Jobs')
        return ('Failure submitting readfibermap Jobs')
    if queue1 is None:
        splog.info('No New Fibermaps Read')
        return (None)
    queue1.monitor_job(pause=float(config.pipe['monitor.pause']), 
                       jobname='slurm_readfibermap')
    return (None)
    
def build_traceflats(mjd, obs, pause=300, 
                     skip_plan=False,
                     allemail=False, **kwrds):
                     
    mjds = ','.join(np.atleast_1d(mjd).astype(str).tolist())

    if config.pipe['Stage.run_spTrace']:
        send_email('build_traceflats '+config.pipe['general.RUN2D'] +' MJD='+mjds +' OBS='+','.join(obs),
                    ptt.join(daily_dir, 'etc','emails'), None,
                    allemail=allemail)
        
    slurm_spTrace.setup_run(alloc=config.queue.get('alloc'),
                            partition=config.queue.get('partition'),
                            nbundle=config.queue.get('nbundle'), mjd=mjd, 
                            nodes=1)
    status = 'Pass'
    attachments = []
    for ob in obs:
        config.pipe['fmjdselect.trace_all_mjds'] = True
        queue1, logfile, errfile = slurm_spTrace.build(mjd, ob)
        if queue1 is None:
            continue
        attachments.extend([logfile,errfile])
        queue1.monitor_job(pause=pause, jobname='slurm_spTrace')
        for mj in np.atleast_1d(mjd):
            for ob in obs:
                plan = read_table_yanny(ptt.join(config.pipe['general.BOSS_SPECTRO_REDUX'],
                                                 config.pipe['general.RUN2D'],'trace',
                                                str(mj),f'spPlanTrace-{mj}_{ob.upper()}.par'),'SPEXP')
                plan.convert_bytestring_to_unicode()
                arcs = plan['flavor'] == 'arc'
                for ff in plan[arcs]['name'].data.flatten():
                    ff = ff.replace('.fit','.fits').replace('sdR','spTraceTab')
                    if not ptt.exists(ptt.join(config.pipe['general.BOSS_SPECTRO_REDUX'],
                                               config.pipe['general.RUN2D'],'trace',str(mj),ff)):
                        status = 'Fail'
    return status, attachments

def build_run(logdir, mj, today, plates=False, traceflat=False):  
    flags = ''
    flags1d = ''
    attachments = None
    obs = config.pipe['fmjdselect.obs']
    epoch = config.pipe['fmjdselect.epoch']
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


    if config.pipe['Stage.run_plan']:
        lco = True if obs[0].upper() == 'LCO' else False
        try:
            spPlan_clobber = config.pipe['Clobber.clobber_plan']
            splog.info('Creating spPlan Files')
            splog.pause_file()
            args = dict(topdir=config.pipe['general.BOSS_SPECTRO_REDUX'], 
                        run2d=config.pipe['general.RUN2D'], 
                        mjd=mj, lco=lco, plates=plates,
                        no_dither=(not config.pipe['fmjdselect.dither']), 
                        returnlist=True, clobber = spPlan_clobber, single_flat = True)
            plans2d = spplan2d(**args)
            
            args = dict(topdir=config.pipe['general.BOSS_SPECTRO_REDUX'], 
                        run2d=config.pipe['general.RUN2D'], 
                        mjd=mj, lco=lco, plates=plates,
                        daily=True, clobber = spPlan_clobber, plans=plans2d)
            spplan1d(**args)
            splog.unpause_file()

        except Exception as e: # work on python 3.x
            splog.error('Failure in building spPlans: '+ str(e))
            if config.pipe['monitor.pipe_monitor']:
                splog.close()
                send_email('Failure '+config.pipe['general.RUN2D'] +' MJD='+mjsub +' OBS='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), mjfile)
                splog.close_file()
            exit
    else:
        plans2d = []
        afc = Field(config.pipe['general.BOSS_SPECTRO_REDUX'], 
                    config.pipe['general.RUN2D'], '*')
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
    no_fibermap = not config.pipe['Stage.run_fibermap']
    if len(plans2d) == 0:
        no_fibermap = True
        traceflat = False
        
    if not no_fibermap:
        splog.info('Building spFibermaps for spplan2ds')

        splog.pause_file()
        error = build_fibermaps(plans2d, mj, obs)
        splog.unpause_file()
        if error is not None:
            splog.close()
            send_email('Failure submitting readfibermap Jobs '+mjsub+' obs='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), [mjfile])
            splog.close_file()
            exit()

    elif no_fibermap:
        splog.info('Skipping pre-Build of spFibermaps for spplan2ds')
        
    args = dict(active=traceflat,
                pause=float(config.pipe['monitor.pause']))        
    if traceflat: #traceflat:
        splog.info('Building TraceFlats for mjd')
        splog.pause_file()
        status, spTatt = build_traceflats(mj, obs, **args)
        splog.unpause_file()
        if status == 'Fail':
            splog.error('Failure in building spTraceFlats and spTraceArcs')
            if config.pipe['monitor.pipe_monitor']:
                attachments = []
                for f in spTatt:
                    if f is None:
                        continue
                    if ptt.exists(f):
                        attachments.append(f)
                if len(attachments) == 0:
                    attachments = None
                splog.close()
                send_email('spTrace Failure '+config.pipe['general.RUN2D'] +' MJD='+mjsub +' OBS='+','.join(obs),
                            ptt.join(daily_dir, 'etc','emails'), attachments)
                splog.close_file()
            exit()
    else:
        splog.info('Skipping Building of TraceFlats for mjd')
        spTatt = [None]
    
    es = '' if not config.pipe['fmjdselect.epoch'] else ' --epoch'
    splog.emailer()


    options = {'run2d': config.pipe['general.RUN2D'],
                'run1d': config.pipe['general.RUN1D'],
                'topdir':config.pipe['general.BOSS_SPECTRO_REDUX'],
                'epoch': config.pipe['fmjdselect.epoch'],
                'mjd': np.atleast_1d(mj).tolist(),
                'obs': obs,
                'MWM_fluxer': config.pipe['reduce.MWM_fluxer'],
                'map3d': config.pipe['reduce.map3d'],
                'no_reject': config.pipe['combine.no_reject'],
                'no_merge_spall': not config.pipe['Stage.run_Sumamrymerge'],
                'a2t': config.pipe['reduce.force_arc2trace'],
                'email': config.pipe['email.sendemail'],
                'allemail': config.pipe['email.allemail'],
                'daily'          : True,
                'clobber'        : config.pipe['Clobber.clobber_pipe'],
                'saveraw': config.pipe['reduce.saveraw'],
                'debug': config.pipe['reduce.debug'],
                'no_db': config.pipe['fibermap.no_db'],
                'no_healpix'     : not config.pipe['Stage.run_healpix'],
                'nodist'         : config.pipe['combine.nodist'],
                }


    pd_ops = pd.Series(options)
    pd_ops = pd.concat([pd.Series(config.queue.to_dict(None))])

    rootlogger.console.setLevel(logging.CRITICAL + 1)
    rootlogger.info('\n'+pd_ops.to_string())
    rootlogger.console.setLevel(logging.DEBUG)
    splog.info('Running uubatchpbs --run2d '+config.pipe['general.RUN2D']+' --obs '+obs[0]+
                     ' --sdssv'+' --email'+
                     ' --topdir '+config.pipe['general.BOSS_SPECTRO_REDUX']+
                     ' --run1d '+config.pipe['general.RUN1D']+ es +
                     f' --sc {config.queue._queue_config}'+
                     ' --mjd '+' '.join(np.asarray(mj).astype(str).tolist()))
    splog.info('')
    
    queue1, redux = uubatchpbs(True)
    
    if config.pipe['monitor.pipe_monitor'] and not config.queue.get('no_write'):
        jobname = f"{config.pipe['general.RUN2D']} MJD={','.join(np.asarray(mj).astype(str).tolist())} OBS={','.join(obs)}"
        if not config.queue.get('nosubmit'):
            subj = queue1.monitor_job(pause=float(config.pipe['monitor.pause']), 
                                      jobname = 'uubatch '+jobname, return_status = True)
        else:
            subj = f'{jobname} not submitted at {datetime.datetime.today().ctime()}'
        if ("not submitted" in subj) or ("Failure" in subj):
            build_summary = False
        build_summary = config.pipe['Stage.run_Summarymerge']
        if build_summary:
            attach_summ, subj_sum = dailysummary( config.pipe['general.module'], mj,
                                                  epoch = config.pipe['fmjdselect.epoch'], 
                                                  pause=float(config.pipe['monitor.pause']), obs=obs,
                                                  jobname = jobname, no_submit=config.queue.get('nosubmit'))
            if attach_summ is None:
                attach_summ = [None]
        else:
            attach_summ = [None]
            subj_sum = None
    else:
        attach_summ = [None]
        subj_sum = None
    splog.close()
    
    if config.pipe['monitor.pipe_monitor'] and not config.queue.get('no_write'):
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
                            topdir=config.pipe['general.BOSS_SPECTRO_REDUX'], 
                            run2d=config.pipe['general.RUN2D'], 
                            run1d=config.pipe['general.RUN1D'], redux = redux)
            
    splog.close_file()
    return
    

def uurundaily():
    if not hasslurm:
        raise(Exception('No slurm package'))

    config.pipe['general.RUN2D'] = load_env('RUN2D')
    config.pipe['general.RUN2D'] = load_env('RUN1D')
    config.pipe['general.BOSS_SPECTRO_REDUX'] = load_env('BOSS_SPECTRO_REDUX')

    obs = config.pipe['fmjdselect.obs']

    if obs[0].lower() != 'apo':
        boss_spectro_data = load_env('BOSS_SPECTRO_DATA_S')
    else:
        boss_spectro_data = load_env('BOSS_SPECTRO_DATA_N')


    if not config.pipe['fmjdselect.epoch']:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')
    else:
        nextmjd_file = ptt.join(daily_dir,'etc','nextmjd_epoch.par')
        flag_file = ptt.join(daily_dir,'etc','completemjd.par')

    today = datetime.datetime.today().strftime("%m%d%Y")
    logdir = ptt.join(daily_dir, "logs", obs[0].upper(),config.pipe['general.RUN2D'])


    makedirs(logdir,exist_ok=True)

    logf = 'uurundaily-'+today+'.log' if not config.pipe['fmjdselect.epoch'] else 'uurunepoch-'+today+'.log'
    rootlogger.open(logfile = ptt.join(logdir, logf), append=True)

    
    rootlogger.debug('========================================')
    rootlogger.debug('Starting at '+datetime.datetime.today().ctime())


    if obs[0].lower() != 'apo': 
        boss_spectro_data = boss_spectro_data.replace('apo', obs[0].lower())
   
    mjd = config.pipe['fmjdselect.mjd']
    if mjd is not None:
        manual=True
    else:
        mjd = get_MJD(boss_spectro_data, config.pipe['general.module'], 
                      obs[0].upper(), config.pipe['general.RUN2D'],
                      nextmjd_file = nextmjd_file, flag_file = flag_file,
                      epoch=config.pipe['fmjdselect.epoch'])
        manual=False
    if len(mjd) > 0:
        if manual is False:
            increment_nextmjd(config.pipe['general.module'], 
                              obs[0].upper(), max(mjd)+1,
                              nextmjd_file = nextmjd_file)

        rootlogger.info('')

        if config.pipe['fmjdselect.batch_mjd'] is True:
            mjd = np.asarray(mjd)
            plate_mjds = mjd[np.where(mjd <  59540)[0]]
            fps_mjds   = mjd[np.where(mjd >= 59540)[0]]

            args = (logdir, plate_mjds.tolist(), today)
            if len(plate_mjds) >0:
                build_run(*args, plates = True, traceflat=False)#**kwrds)

            if len(fps_mjds) > 0:
                args = (logdir, fps_mjds.tolist(), today)
                build_run(*args, plates=False, traceflat=True)#**kwrds)

        else:
            for mj in mjd:
                if mj < 59540:
                    plates = True
                    rtf= False
                else:
                    plates = False
                    rtf = config.pipe['Stage.run_spTrace']
                args = (logdir, [mj], today)
                build_run(*args, plates=plates, traceflat = rtf)#**kwrds)

        if (not manual) and config.pipe['monitor.pipe_monitor']:
            flag_complete(config.pipe['general.module'], mjd, obs[0].upper(), flag_file = flag_file)
    rootlogger.debug('Completed at '+datetime.datetime.today().ctime())


def parseNumList(string):
    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise argparse.ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start), int(end)+1))


class SkipPlanAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values is None:
            setattr(namespace, self.dest, True)
        elif values == 'all':
            setattr(namespace, self.dest, True)
        else:
            setattr(namespace, self.dest, values)
