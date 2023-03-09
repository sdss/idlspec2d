#!/usr/bin/env python3
import argparse
from os import getenv, makedirs, popen, chdir, getcwd, popen
import os.path as ptt
from uubatchpbs import uubatchpbs, Formatter
from load_module import load_module
from pydl.pydlutils.yanny import yanny, write_table_yanny
import numpy as np
from astropy.table import Table
import logging
import datetime
import astropy.time
import sys
from glob import glob
import re

jdate = str(int(float(astropy.time.Time(datetime.datetime.utcnow()).jd) - 2400000.5))

def email(subj='', message='', email_file = ptt.join(getenv('HOME'),'daily','etc','emails')):
    if ptt.exists(email_file):
        for address in open(email_file).read().splitlines():
            if len(address) == 0: 
                continue
            stream = popen('echo "'+message+'"| mail -s "'+subj+'" '+address)
            output = stream.read()
    else:
        logging.warning('no file at '+email_file)
        

def printAndRun(log, cmd, idlspec2d_dir):
    log.info('Running '+cmd)
    cmd.replace('spplan', idlspec2d_dir+'/pro/spec2d/spplan')
    stream = popen(cmd)
    log.info(stream.read())
    log.info('')

def load_env(key):
    val = getenv(key)
    if val is None:
        print('ERROR: '+key+' is not set')
        exit()
    return(val)

def read_module(mod):
    module = load_module()
    module('purge')
    module('load', mod)
    run2d = load_env('RUN2D')
    run1d = load_env('RUN1D')
    boss_spectro_redux = load_env('BOSS_SPECTRO_REDUX')
    boss_spectro_data_N =load_env('BOSS_SPECTRO_DATA_N')
    idlspec2d_dir =load_env('IDLSPEC2D_DIR')
    return(run2d, run1d, boss_spectro_redux, boss_spectro_data_N,idlspec2d_dir)

def get_nextmjd(mod, obs, nextmjd_file = ptt.join(getenv('HOME'),'daily','etc','nextmjd.par')):
    nextmjds = yanny(nextmjd_file)

    obss  = np.char.upper(nextmjds["NEXTMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) == 0:
        nextmjd = jdate
    else:
        nextmjd = nextmjds["NEXTMJD"]['mjd'][indx][0]
    return(nextmjd)

def increment_nextmjd(logger, mod, obs, nextmjd, nextmjd_file = ptt.join(getenv('HOME'),'daily','etc','nextmjd.par')):
    nextmjds = yanny(nextmjd_file)
    obss  = np.char.upper(nextmjds["NEXTMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) != 0: nextmjds["NEXTMJD"]['mjd'][indx] = nextmjd
    tab_nextmjds = Table(nextmjds["NEXTMJD"])
    if len(indx) == 0: tab_nextmjds.add_row([mod, nextmjd, obs])
    write_table_yanny(tab_nextmjds, nextmjd_file, tablename = "NEXTMJD", overwrite = True)
    logger.info("Next MJD to wait for will be "+str(nextmjd))

def get_MJD(logger, boss_spectro_data, mod, obs, run2d, nextmjd_file = ptt.join(getenv('HOME'),'daily','etc','nextmjd.par')):
    nextmjd = get_nextmjd(mod, obs, nextmjd_file = nextmjd_file)
    logger.info("Looking for MJDs of or after "+str(nextmjd))
    path = ptt.join(boss_spectro_data, '?????')
    def get_key(fp):
        if not ptt.isdir(fp): return(0)
        filename = ptt.basename(fp)
        int_part = filename.split()[0]
        return(int(int_part))
    files = sorted(glob(path),key=get_key)
    lastmjd = int(ptt.basename(files[-1]))
    mjd = []
    while lastmjd >= int(nextmjd):
        if ptt.isdir(path.replace('?????', str(lastmjd))):
            mjd.append(lastmjd)
        else:
            email(subj = 'skipping '+str(lastmjd)+' for '+mod+' obs='+obs)
        lastmjd = lastmjd - 1
    if len(mjd) == 0:
        logger.info('MJD '+str(nextmjd)+' for run2d='+run2d+' OBS='+obs+' is not here yet')
    else:
        logger.info('MJDs for run2d='+run2d+' OBS='+obs+ ' transfered: '+str(nextmjd)) 
    return mjd

def build_run(skip_plan, logdir, obs, mj, run2d, run1d, idlspec2d_dir, options, topdir, today, plates = False):
    flags = ''
    if plates is True:
        flags = flags + ', /plate_s'
    if obs[0].upper() == 'LCO':
        flags = flags + ', /lco'

    if len(mj) == 1:
        mjfilelog = logging.FileHandler(ptt.join(logdir, str(mj[0])+'.log'))
    else:
        mjfilelog = logging.FileHandler(ptt.join(logdir, str(mj[0])+'-'+str(mj[-1])+'.log'))
    mjfilelog.setLevel(logging.DEBUG)
    mjfilelog.setFormatter(Formatter())
    mjconsole = logging.StreamHandler()
    mjconsole.setLevel(logging.DEBUG)
    mjconsole.setFormatter(Formatter())

    rootfilelog = logging.FileHandler(ptt.join(logdir, 'uurundaily-'+today+'.log'))
    rootfilelog.setLevel(logging.DEBUG)
    rootfilelog.setFormatter(Formatter())
    logger = logging.getLogger(str(mj))
    logger.addHandler(rootfilelog)
    logger.addHandler(mjfilelog)
    logger.addHandler(mjconsole)
    logger.setLevel(logging.DEBUG)
    if not skip_plan:
        for mjd in mj:
            printAndRun(logger, "idl -e 'spplan2d"+flags+", MJD="+str(mjd)+"'",idlspec2d_dir)
            printAndRun(logger, "idl -e 'spplan1d"+flags+", MJD="+str(mjd)+"'",idlspec2d_dir)
    else:
        logger.info('Using old spplan files')
    logger.info('Running uubatchpbs.py --run2d '+run2d+' --obs '+obs[0]+' --sdssv_fast --email'+
                     ' --topdir '+topdir+ ' --run1d '+run1d+
                     ' --mjd '+' '.join(np.asarray(mj).astype(str).tolist()))
    logger.info('')
    uubatchpbs(**options, obs=obs, run2d = run2d, run1d = run1d, topdir = topdir, mjd=mj, logger=logger)
    
    logger.removeHandler(mjconsole)
    logger.removeHandler(rootfilelog)
    logger.removeHandler(mjfilelog)
    mjfilelog.close()
    mjconsole.close()
    rootfilelog.close()


def uurundaily(module, obs, mjd = None, clobber=False, fast = False, saveraw=False, skip_plan=False, 
              nosubmit=False, noslurm=False, batch=False, debug=False, nodb=False):
    run2d, run1d, topdir, boss_spectro_data, idlspec2d_dir = read_module(module)
    nextmjd_file = ptt.join(getenv('HOME'),'daily','etc','nextmjd.par')
    today = datetime.datetime.today().strftime("%m%d%Y")
    logdir = ptt.join(getenv("HOME"), "daily", "logs", obs[0].upper(),run2d.upper())


    makedirs(logdir,exist_ok=True)
    rootlogger = logging.getLogger('root')

    filelog = logging.FileHandler(ptt.join(logdir, 'uurundaily-'+today+'.log'))
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
        mjd = get_MJD(rootlogger, boss_spectro_data, module, obs[0].upper(), run2d, nextmjd_file = nextmjd_file)
        manual=False
    if len(mjd) > 0:
        if manual is False:
            increment_nextmjd(rootlogger, module, obs[0].upper(), max(mjd)+1, nextmjd_file = nextmjd_file)

        options = {'MWM_fluxer'     : True,
                   'no_reject'      : True,
                   'no_merge_spall' : True,
                   'walltime'       : '40:00:00',
                   'include_bad'    : True,
                   'xyfit'          : True,
                   'loaddesi'       : True,
                   'shared'         : True,
                   'mem_per_cpu'    : '7500',
                   'fast'           : fast,
                   'email'          : True,
                   'nosubmit'       : nosubmit,
                   'daily'          : True,
                   'clobber'        : clobber,
                   'saveraw'        : saveraw,
                   'debug'          : debug,
                   'no_db'          : nodb,
                   'no_write'       : noslurm,
                   }
        rootlogger.info('')
        if batch is True:
            mjd = np.asarray(mjd)
            plate_mjds = mjd[np.where(mjd <  59540)[0]]
            if len(plate_mjds) >0:
                build_run(skip_plan, logdir, obs, plate_mjds.tolist(), run2d, run1d, idlspec2d_dir, options, topdir, today, plates = True)
            fps_mjds   = mjd[np.where(mjd >= 59540)[0]]
            if len(fps_mjds) > 0:
                build_run(skip_plan, logdir, obs, fps_mjds.tolist(), run2d, run1d, idlspec2d_dir, options, topdir, today, plates = False)
        else:
            for mj in mjd:
                if mj < 59540:
                    plates = True
                else: 
                    plates = False
                build_run(skip_plan, logdir, obs, [mj], run2d, run1d, idlspec2d_dir, options, topdir, today, plates = plates)


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
            nosubmit=args.nosubmit, batch=args.batch, noslurm=args.noslurm, debug=args.debug, nodb= args.nodb)
