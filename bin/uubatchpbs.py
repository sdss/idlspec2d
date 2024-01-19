#!/usr/bin/env python3

from slurm import queue
import argparse
from os import getenv, makedirs, popen, chdir, getcwd, popen
import os.path as ptt
from glob import glob
from pydl.pydlutils.yanny import yanny, read_table_yanny
import sys
import numpy as np
import io
import datetime
import astropy.time
import pandas as pd
import logging
from dailylogger import *
import time

if getenv('SLURM_VERS') == 'notchpeak': 
    share = True
else: 
    share = False
jdate = int(float(astropy.time.Time(datetime.datetime.utcnow()).jd)-2400000.5)

def get_dirs(basedir, val=None, start=None, end=None, format = '??????', numeric=True):
    dirs = glob(ptt.join(basedir, format))
    valid_dirs = []
    if val is not None:
        for i,v in enumerate(val):
            val[i] = str(v).zfill(len(format))

    for fd in dirs:
        if numeric is True:
            if not ptt.basename(fd).isnumeric():
                continue
        if val is not None:
            if ptt.basename(fd) not in val:
                continue
        if start is not None:
            if int(ptt.basename(fd)) < int(start): 
                continue
        if end is not None:
            if int(ptt.basename(fd)) > int(end):
                continue
        valid_dirs.append(fd)

    return(valid_dirs)

def mjd_match(thismjd, mjd=None, mjdstart=None, mjdend=None):
    if mjd is not None:
        if int(thismjd) not in np.atleast_1d(np.asarray(mjd)).astype(int).tolist():
            return(False)
    if mjdstart is not None:
        if int(thismjd) < int(mjdstart):
            return(False)
    if mjdend is not None:
        if int(thismjd) > int(mjdend):
            return(False)
    return(True)


def build_cmd(topdir=None,run2d=None,run1d=None,idlutils_1d=None,
        no_reject=False, MWM_fluxer=False, map3d='bayestar15',
        noxcsao=False, skip_specprimary=False,
        no_merge_spall=False, skip2d=False, onestep_coadd=False,
        fibermap_clobber=False, lco=False, plan2d=None, plancomb=None,
        fieldmjd=None, post_idl=False, only1d=False, daily=False, module="",
        custom=None, allsky=False, epoch=False, saveraw=False, debug=False,
        sdss_access_release = None, sdss_access_remote = False,
        no_db=False, fast_no_db=False,no_healpix=False, dr19=False,
        custom_coadd_only=False, custom_1dpost=False, redux=None, **kwargs):

    field = fieldmjd.split('-')[-2]
    mjd = fieldmjd.split('-')[-1]

    cmd = []
    cmd.append('# Auto-generated batch file '+datetime.datetime.now().strftime("%c"))
    if daily:
        cmd.append(f"module purge ; module load {module}")
        cmd.append('')

    if epoch is False:
        if custom is None:
            cmd.append('cd '+ptt.join(topdir,run2d,field))
        else:
            if allsky is False:
                cmd.append('cd '+ptt.join(topdir,run2d,custom,field))
                fieldmjd = custom+'-'+mjd
            else:
                cmd.append('cd '+ptt.join(topdir,run2d,custom))
                fieldmjd = custom+'-'+mjd
    else:
        if custom is None:
            cmd.append('cd '+ptt.join(topdir,run2d,field,'epoch'))
        else:
            if allsky is False:
                cmd.append('cd '+ptt.join(topdir,run2d,custom,field,'epoch'))
                fieldmjd = custom+'-'+mjd
            else:
                cmd.append('cd '+ptt.join(topdir,run2d,custom,'epoch'))
                fieldmjd = custom+'-'+mjd
    if not allsky:
        if lco:
            cmd.append('export BOSS_SPECTRO_DATA=$BOSS_SPECTRO_DATA_S')
            cmd.append('export GCAM_DATA=$GCAM_DATA_S')
        else:
            cmd.append('export BOSS_SPECTRO_DATA=$BOSS_SPECTRO_DATA_N')
            cmd.append('export GCAM_DATA=$GCAM_DATA_N')
    
    cmd.append('')
    cmd.append("#- Echo commands to make debugging easier")
    cmd.append("set -o verbose")

    spreduce2d_keys=''
    rm_combine_keys='/xyfit, /loaddesi,'
    spec1d_keys=''
    xcsao_keys=''
    
    spcalib_keywords = ''
    if post_idl:
        reformat_keywords = ''
        conflist_keywords = ' /create,'
        fieldmerge_keywords = ' /include_bad,'
        image_keywords = ''
    else:
        spSpecRef_key = "" # " --lsdr10"
        flist_key = ""
        fmerge_key = " --include_bad"

    if epoch:
        spSpecRef_key    = spSpecRef_key    + " --epoch"
        flist_key        = flist_key        + " --epoch"
        fmerge_key       = fmerge_key       + " --epoch"
        rm_combine_keys  = rm_combine_keys  + " /epoch,"
        spec1d_keys      = spec1d_keys      + " /epoch,"
        xcsao_keys       = xcsao_keys       + " --epoch"
        spcalib_keywords = spcalib_keywords + " /epoch,"
        
    fmap_keywords = ''
    if fibermap_clobber:
        fmap_keywords = fmap_keywords + " --clobber"
    if no_db:
        fmap_keywords = fmap_keywords + " --no_db"
        if fast_no_db:
            fmap_keywords = fmap_keywords + " --fast"
    if sdss_access_release is not None:
        fmap_keywords = fmap_keywords + " --release "+sdss_access_release
    if sdss_access_remote is not None:
        fmap_keywords = fmap_keywords + " --remote"

    if MWM_fluxer:
        spreduce2d_keys = spreduce2d_keys + ' /MWM_fluxer,'
        rm_combine_keys = rm_combine_keys + ' /MWM_fluxer,'
        spreduce2d_keys = spreduce2d_keys + f' map3d={map3d},'
    if saveraw:
        spreduce2d_keys = spreduce2d_keys + ' /saveraw,'
    if debug:
        spreduce2d_keys = spreduce2d_keys + ' /debug,'
    if no_reject:
        rm_combine_keys = rm_combine_keys + ' /no_reject,'
    if onestep_coadd:
        rm_combine_keys = rm_combine_keys + ' /onestep_coadd,'

    
    if post_idl:
        if not noxcsao:
            reformat_keywords = reformat_keywords+' /XCSAO,'
            fieldmerge_keywords = fieldmerge_keywords+' /XCSAO,'
        if skip_specprimary:
            fieldmerge_keywords = fieldmerge_keywords+' /skip_specprimary, '
    else:
        if not noxcsao:
            fmerge_key = fmerge_key +" --XCSAO"
        if skip_specprimary:
            fmerge_key = fmerge_key +" --skip_specprimary"


    if custom is not None:
        customkey = ' custom='+custom+','
        pycustomkey = ' --custom '+custom
        if allsky: 
            customkey = customkey + ' /allsky,'
            pycustomkey = ' --allsky'
            spcalib_keywords    = spcalib_keywords + ' /allsky,'

        rm_combine_keys     = rm_combine_keys     + customkey
        spec1d_keys         = spec1d_keys         + customkey
        #xcsao_keys          = xcsao_keys          + pycustomkey
        if post_idl:
            reformat_keywords   = reformat_keywords   + customkey
            conflist_keywords   = conflist_keywords   + customkey
            fieldmerge_keywords = fieldmerge_keywords + customkey
            image_keywords      = image_keywords      + customkey
        else:
            fmerge_key = fmerge_key + pycustomkey
            flist_key = flist_key + pycustomkey
            spSpecRef_key = spSpecRef_key + pycustomkey

    if epoch:
        skip2d = True
    if only1d:
        skip2d = True
    cmd.append('')
    cmd.append('#- The real work')
    

    if not skip2d:
        for plan in plan2d:
            plan = plan.strip("'")
            if not dr19:
                cmd.append(f"readfibermaps.py --spplan2d {plan}")
            else:
                cmd.append(f"readfibermaps.py --spplan2d {plan} --dr19")
            cmd.append('touch '+plan.replace('.par', '.started').replace('spPlan2d','spec2d'))
            cmd.append("echo 'spreduce2d,"+spreduce2d_keys+' "'+plan+'"'+"' | idl")
            cmd.append('touch '+plan.replace('.par', '.done').replace('spPlan2d','spec2d'))
    if not custom_1dpost:
        if custom is not None:
            for plan in plan2d:
                plan = plan.strip("'")
                cmd.append('touch '+plan.replace('.par', '.started').replace('spPlanCustom','specombine'))
                cmd.append("echo 'spspec_target_merge, "+' "'+plan+'"'+"' | idl")
                cmd.append('touch '+plan.replace('.par', '.done').replace('spPlanCustom','specombine'))
            if custom_coadd_only:
                return(cmd)
    
    if not only1d:
        if epoch:
            cmd.append('touch spPlancombepoch-'+fieldmjd+'.started')
            cmd.append("echo 'rm_combine_script, "+'"spPlancombepoch-'+fieldmjd+'.par",'+rm_combine_keys+' run2d="'+run2d+'"'+"' |idl")
            cmd.append('touch spPlancombepoch-'+fieldmjd+'.done')
        elif custom is None:
            cmd.append('touch specombine-'+fieldmjd+'.started')
            cmd.append("echo 'rm_combine_script, "+'"spPlancomb-'+fieldmjd+'.par",'+rm_combine_keys+' run2d="'+run2d+'"'+"' |idl")
            cmd.append('touch specombine-'+fieldmjd+'.done')


    if run1d is None: 
        run1d=run2d
    if run2d != run1d: 
        cmd.append('module switch idlspec2d idlspec2d/'+run1d)
    if idlutils_1d is not None:
        cmd.append('module switch idlutils idlutils/'+idlutils_1d)

    if allsky is False:
        cmd.append('touch spec1d-'+fieldmjd+'.started')
        cmd.append("echo 'spreduce1d_empca, "+'"spField-'+fieldmjd+'.fits",'+spec1d_keys+' run1d="'+run1d+'"'+"' |idl")
        if not noxcsao:
            cmd.append('run_PyXCSAO.py spField-'+fieldmjd+'.fits'+xcsao_keys+' --run1d "'+run1d+'"')
        cmd.append('touch spec1d-'+fieldmjd+'.done')
    else:
        global plancomb_last
        global EPOCH_COMBINE
        if 'plancomb_last' not in globals():
            plancomb_last = ''
        if plancomb != plancomb_last:
            cplan = read_table_yanny(plancomb, 'COADDPLAN')
            EPOCH_COMBINE = np.sort(np.unique(cplan['EPOCH_COMBINE']))
            plancomb_last = plancomb
        for mjec in EPOCH_COMBINE:
            if custom_1dpost:
                test_mjec = redux.split('_')[-1]
                if str(mjec) != test_mjec:
                    continue
            fmjd = custom+'-'+str(mjec)
            cmd.append('touch spec1d-'+fmjd+'.started')
            cmd.append("echo 'spreduce1d_empca, "+'"spFullsky-'+fmjd+'.fits",'+spec1d_keys+' run1d="'+run1d+'"'+"' |idl")
            if not noxcsao:
                cmd.append('run_PyXCSAO.py spFullsky-'+fmjd+'.fits'+xcsao_keys+' --run1d "'+run1d+'"')
            cmd.append('touch spec1d-'+fmjd+'.done')


    if not only1d:
        if post_idl:
            cmd.append('')
            cmd.append('#- Make final spectra files')
            cmd.append("echo 'reformat_spec, "+'"spField-'+fieldmjd+'.fits",'+reformat_keywords+' run1d="'+run1d+'"'+"' |idl")
            cmd.append("echo 'conflist,"+conflist_keywords[:-1]+"' |idl")
            cmd.append("echo 'fieldmerge, field="+field+", mjd="+mjd+","+fieldmerge_keywords[:-1]+"' |idl")
            cmd.append("echo 'reformat_spec, "+'"spField-'+fieldmjd+'.fits", /lite,'+reformat_keywords+' run1d="'+run1d+'"'+"' |idl")

            cmd.append('')
            cmd.append('#- Make Spectro-Photometric Calibration QA Figures')
            cmd.append("echo 'spcalib_qa, fieldid="+field+", mjd="+mjd+', '+spcalib_keywords+' run2d="'+run2d+'"'+"' |idl")

            cmd.append('')
            cmd.append('#- Make pretty pictures')
            cmd.append("echo 'plate_spec_image, "+field+", mjd="+mjd+', '+image_keywords+' run1d="'+run1d+'", run2d="'+run2d+'"'+"' |idl")

            if not no_merge_spall:
                cmd.append('#- update spAll file')
                cmd.append("echo 'fieldmerge, "+fieldmerge_keywords[:-1]+"' |idl")

        else:
            if allsky is False:
                fmjd = [fieldmjd]
            else:
                fmjd = []
                for mjec in EPOCH_COMBINE:
                    if custom_1dpost:
                        test_mjec = redux.split('_')[-1]
                        if str(mjec) != test_mjec:
                            continue
                    fmjd.append( custom+'-'+str(mjec))
            for fm in fmjd:
                field, mjd = fm.split('-')
                if allsky is False:
                    cmd.append('')
                    cmd.append('#- Make Field List')
                    cmd.append(f"fieldlist.py --create --run1d {run1d} --run2d {run2d} {flist_key} --logfile fieldlist-{field}-{mjd}.log")
                
                    cmd.append('')
                    cmd.append('#- Make Field spAll file')
                    cmd.append(f"fieldmerge.py {fmerge_key} --clobber --run2d {run2d} --field {field} --mjd {mjd}")
                else:
                    cmd.append('')
                    cmd.append('#- Make Field spAll file')
                    cmd.append(f"fieldmerge.py {fmerge_key} --clobber --run2d {run2d} --custom {field} --run1d {run1d} --mjd {mjd}")
                cmd.append('')
                cmd.append('#- Make final spectra files')
                cmd.append(f"spSpec_reformat.py --field {field} --mjd {mjd} -p --run1d {run1d} --run2d {run2d} {spSpecRef_key}")

                if allsky is False:
                    cmd.append('')
                    cmd.append('#- Make Spectro-Photometric Calibration QA Figures')
                    cmd.append("echo 'spcalib_qa, fieldid="+field+", mjd="+mjd+','+spcalib_keywords+' run2d="'+run2d+'"'+"' |idl")

            if not no_merge_spall:
                cmd.append('')
                cmd.append('#- Making Summary Files')
                if allsky is False:
                    cmd.append(f"fieldlist.py --create --run1d {run1d} --run2d {run2d} {flist_key}")
                    cmd.append(f"fieldmerge.py --lite {fmerge_key} --run2d {run2d} --remerge_fmjd {field}-{mjd}")
                else:
                    cmd.append(f"fieldmerge.py --lite {fmerge_key} --custom {field} --run1d {run1d} --run2d {run2d} --remerge_fmjd {field}-{mjd}")
                    
        if custom is None:
            if not no_healpix:
                cmd.append('#- Make the healpix links')
                cmd.append('sas_mwm_healpix --spectro boss --mjd '+mjd+' --telescope apo25m --drpver '+run2d+' -v')

    return(cmd)

def uubatchpbs(obs = ['apo', 'lco'], topdir = getenv('BOSS_SPECTRO_REDUX'),
               run1d = getenv('RUN1D'), run2d = getenv('RUN2D'), idlutils_1d = None,
               no_reject = False, MWM_fluxer = True, map3d='bayestar15',
               noxcsao = False, skip_specprimary = False,
               no_merge_spall = False, skip2d = False, onestep_coadd = False,
               fibermap_clobber = False,
               include_bad = False, xyfit = False, loaddesi = False, only1d = False,
               fieldids = None, fieldstart = None, fieldend = None,
               mjd = None, mjdstart = None, mjdend = None, saveraw=False,
               no_write = False, kingspeak = False, shared = share, fast = False,
               mem_per_cpu = getenv('SLURM_MEM_PER_CPU'), walltime = '336:00:00',
               nodes = None, ppn = None, nosubmit = False, daily=False, module="",
               debug=False, no_db=False, dr19=False,
               clobber = False, custom= None, allsky = False, epoch = False, no_healpix=False,
               email = False, logger=None, fast_no_db=False,
               sdss_access_remote = False, sdss_access_release=None,
               custom_coadd_only=False, custom_1dpost=False):


    elog = emailLogger()
    emaillog = elog.log_handler
    emaillog.setLevel(logging.DEBUG)
    emaillog.setFormatter(Formatter())


    if logger is None:
        logger = logging.getLogger()
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        console.setFormatter(Formatter())
        logger.addHandler(console)
        logger.addHandler(emaillog)
        logger.setLevel(logging.DEBUG)
    else:
        logger.addHandler(emaillog)
        logger.setLevel(logging.DEBUG)
        
    if kingspeak:
        shared=False
    
    cmdinputs = locals()
    fullinputs = cmdinputs.copy()
    cmdinputs.pop('topdir')
    cmdinputs.pop('no_write')
    cmdinputs.pop('logger')
    cmdinputs.pop('daily')
    try:
        cmdinputs.pop('console')
    except:
        pass
    cmdinputs.pop('elog')
    cmdinputs.pop('emaillog')
    if cmdinputs['custom'] is None:
        cmdinputs.pop('allsky')
    for key in list(cmdinputs.keys()):
        if cmdinputs[key] is None:
            cmdinputs.pop(key)
    error = False
    
    if not ptt.isdir(topdir):
        logger.warning('topdir (BOSS_SPECTRO_REDUX) is invalid')
        error = True
    else:
        logger.info('topdir   '+topdir)
    logger.info(pd.Series(cmdinputs).to_string())
    
    topdir2d = ptt.join(topdir, run2d)
    if custom is not None:
        topdir2d = ptt.join(topdir2d, custom)
    if not allsky:
        fielddirs = get_dirs(topdir2d, val = fieldids, start=fieldstart, end=fieldend)
        if len(fielddirs) == 0:
            logger.warning('No Directories Found')
            error = True
    else:
        fielddirs = ['']

    if error:
        logger.removeHandler(emaillog)
        emaillog.close()
        return()

    redux_list = []
    skipped = 0
   
    if epoch is False:
        if custom is None:
            plan_str = 'spPlancomb-*.par'
        else:
            plan_str = 'spPlanCustom-'+custom+'-*.par'
    else:
        if custom is None:
            plan_str = 'spPlancombepoch-*.par'
        else:
            plan_str = 'spPlancombepoch_'+custom+'-*.par'

    for fielddir in fielddirs:
        if epoch is False:
            if custom is None:
                planfile = glob(ptt.join(topdir2d, fielddir, plan_str))
            else:
                planfile = glob(ptt.join(topdir2d, fielddir, plan_str))
        else:
            if custom is None:
                planfile = glob(ptt.join(topdir2d, fielddir, 'epoch', plan_str))
            else:
                planfile = glob(ptt.join(topdir2d, custom, fielddir,'epoch', plan_str))
        ifile = len(planfile)
        for plan in planfile:
            yplan = yanny(plan)
            hdr = yplan.new_dict_from_pairs()
            thismjd = hdr['MJD']
            if 'OBS' in hdr.keys():
                thisobs = hdr['OBS'].lower()
            else:
                thisobs = 'apo'
            if thisobs.lower() not in obs:
                continue
            if not mjd_match(thismjd, mjd=mjd, mjdstart=mjdstart, mjdend=mjdend):
                continue
            #fieldmjd = hdr['fieldid']+'-'+hdr['MJD']
            if thisobs.lower() == 'lco':
                lco = True
            else:
                lco = False

            if epoch is False:
                if custom is None:
                    fieldmjd = hdr['fieldid']+'-'+hdr['MJD']
                    redux = ptt.join(topdir2d, fielddir, 'redux-'+fieldmjd)
                else:
                    if allsky:
                        fieldmjd = custom+'-'+hdr['CreateMJD']
                    else:
                        fieldmjd = custom+'-'+hdr['fieldid']+'-'+hdr['CreateMJD']
                    if not custom_1dpost:
                        redux = ptt.join(topdir2d, fielddir, 'redux_'+fieldmjd)
                    else:
                        redux = []
                        cplan = read_table_yanny(plan, 'COADDPLAN')
                        EPOCH_COMBINE = np.sort(np.unique(cplan['EPOCH_COMBINE']))
                        for mjec in EPOCH_COMBINE:
                            redux.append(ptt.abspath(ptt.join(topdir2d, fielddir,'redux_'+fieldmjd+'_'+str(mjec))))
            else:
                if custom is None:
                    fieldmjd = hdr['fieldid']+'-'+hdr['MJD']
                    redux = ptt.join(topdir2d, fielddir,'epoch','redux-'+fieldmjd)
                    custom_1dpost = False
                else:
                    if allsky:
                        fieldmjd = custom+'-'+hdr['CreateMJD']
                    else:
                        fieldmjd = custom+'-'+hdr['fieldid']+'-'+hdr['CreateMJD']
                    if not custom_1dpost:
                        redux = ptt.join(topdir2d, fielddir,'epoch','redux_'+fieldmjd)
                    else:
                        redux = []
                        cplan = read_table_yanny(plan, 'COADDPLAN')
                        EPOCH_COMBINE = np.sort(np.unique(cplan['EPOCH_COMBINE']))
                        for mjec in EPOCH_COMBINE:
                            redux.append(ptt.abspath(ptt.join(topdir2d, fielddir,'epoch','redux_'+fieldmjd+'_'+str(mjec))))
                            print(redux[-1])
            if not custom_1dpost:
                redux = [ptt.abspath(redux)]
            
            if custom is None:
                plan2d = hdr['planfile2d'].split(' ')
            else:
                plan2d = [ptt.basename(plan)]
        

            for redux1 in redux:
                if (not ptt.exists(redux1)) or (clobber is True):
                    cmd = build_cmd(**fullinputs, plan2d=plan2d, plancomb=plan, fieldmjd=fieldmjd, lco=lco, redux=redux1)
                    if not no_write:
                        with open(redux1,'w') as r:
                            for c in cmd:
                                r.write(c+'\n')
                    redux_list.append(redux1)
                else:
                    skipped += 1
    if len(redux_list) > 25:
        fast=False
    logger.info('')
    logger.info('---------------------------------------------------')
    logger.info('boss_redux: #BOSS Field-MJDs Done  = '+str(skipped))
    logger.info('            #BOSS Field-MJDs To Do = '+str(len(redux_list)))
    logger.info('---------------------------------------------------')
    logger.info('')
    
    if len(redux_list) == 0: 
        logger.removeHandler(emaillog)
        emaillog.close()
        return None


    if kingspeak is False:
        qos = 'sdss'
        alloc = 'sdss-np'
        if fast is True: alloc = alloc+'-fast'
        partition = 'sdss-np'
        max_c = 64
        if ppn is None: ppn= np.min([64, np.max([len(redux_list),2])])
    else:
        qos = 'sdss'
        alloc = 'sdss-kp'
        partition = 'sdss-kp'
        shared=False
        max_c = 16
        if ppn is None: ppn= np.min([16, np.max([len(redux_list),2])])

    if daily is True:
        if nodes is None and len(redux_list) > ppn:
            nodes = int(np.ceil(len(redux_list)/max_c))
        else:
            nodes = 1
    obsstr = '_'.join(obs)
    if custom is None:
        if not epoch:
            label = run2d + '_' + obsstr
        else:
            label = run2d + '_epoch_' + obsstr
    else:
        if not epoch:
            label = run2d + '_' + obsstr + '_'+custom
        else:
            label = run2d + '_epoch_' + obsstr + '_'+custom

    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    if not no_write:
        queue1 = queue(key=None, verbose=True)
        queue1.create(label=label, nodes=str(nodes), ppn=str(ppn), partition = partition,
                      alloc=alloc, shared=shared, walltime=walltime, mem_per_cpu=mem_per_cpu)
    else:
        queue1 = None
    for i in range(len(redux_list)):
        cmd, log, err = make_run_cmd(redux_list[0])
        redux_list.pop(0)
        if not no_write:
            queue1.append(cmd, outfile = log, errfile = err)
        else: 
            logger.info(cmd)
    if not no_write:
        queue1.commit(hard=True, submit=(not nosubmit))
    output = new_stdout.getvalue()
    sys.stdout = old_stdout
    logger.info(output)

    if email is True:
        if daily:
            jdate = mjd[0]
        else:
            jdate = int(float(astropy.time.Time(datetime.datetime.utcnow()).jd)-2400000.5)
        elog.send('UUBATCH '+run2d +' MJD='+str(jdate) +' OBS='+','.join(obs), ptt.join(getenv('HOME'), 'daily', 'etc','emails'), logger)
    logger.removeHandler(emaillog)
    emaillog.close()

    return(queue1)

def make_run_cmd(redux):
    cmd = 'source '+redux
    log = redux+'.o'
    err = redux+'.e'
    return(cmd, log, err)




if __name__ == '__main__' :
    """
    Batch process Spectro-2D and Spectro-1D reductions based upon already-built plan files
    """
    parser = argparse.ArgumentParser(
            prog=ptt.basename(sys.argv[0]),
            description='Build idlspec2d redux and submit to slurm')


    shortgroup = parser.add_argument_group('Short cuts')
    shortgroup.add_argument('--sdssv', action='store_true', help='--mwm, --no_merge_spall, --no_reject')
    shortgroup.add_argument('--sdssv_fast', action='store_true')
    shortgroup.add_argument('--sdssv_noshare', action='store_true')
    shortgroup.add_argument('--apo', action = 'store_true')
    shortgroup.add_argument('--lco', action = 'store_true')
    shortgroup.add_argument('--bay15', action='store_true', help='Set map3d to bayestar15 model')
#    shortgroup.add_argument('--eden23', action='store_true', help='Set map3d to edenhofer2023 model')
    shortgroup.add_argument('--merge3d', action='store_true', help='Set map3d to best 3d model')

    rungroup = parser.add_argument_group('idlspec2d Run options')
    rungroup.add_argument('--obs', help='Observatory {apo,lco}', nargs='*', default=['apo','lco'], type=str.lower)
    rungroup.add_argument('--topdir', type=str, help='Optional override value for the environment variable $BOSS_SPECTRO_REDUX',
                          default = getenv('BOSS_SPECTRO_REDUX'))
    rungroup.add_argument('--run1d', type=str, help='Optional override value for the enviro variable $RUN1D', default=getenv('RUN1D'))
    rungroup.add_argument('--run2d', type=str, help='Optional override value for the enviro variable $RUN2D', default=getenv('RUN2D'))
    rungroup.add_argument('--idlutils_1d', type=str, help='idlutils override version of spec1d', default=None)
    rungroup.add_argument('--no_reject', action='store_true', help='Deactivate Rejection in Coadd')
    rungroup.add_argument('--MWM_fluxer','--mwm', action='store_true', help='')
    rungroup.add_argument('--map3d',type=str.lower,default='bayestar15',
                            help='Name of 3d dustmap to use with MWM_fluxer (default=bayestar15)',
                            choices = ['bayestar15','bay15','merge3d']) #['bayestar15','bay15','edenhofer2023','eden23','merge3d'])
    rungroup.add_argument('--no_healpix','--nohp', action='store_true', help='Turn off copy to healpix')
    rungroup.add_argument('--noxcsao', action='store_true', help='Skip pyXCSAO')
    rungroup.add_argument('--skip_specprimary', action='store_true', help='Skip Calculation of Specprimary')
    rungroup.add_argument('--no_merge_spall', action='store_true', help='Skip building full SpAll File')
    rungroup.add_argument('--skip2d', action='store_true', help='Skip spreduce2d')
    rungroup.add_argument('--only1d', action='store_true', help='run spec1d step only (eg. spreduce1d_empca, XCSAO)')
    rungroup.add_argument('--onestep_coadd', action='store_true', help='Use legacy one step version of coadd')
    rungroup.add_argument('--fibermap_clobber', action='store_true', help='Clobber spfibermap fits file')
    rungroup.add_argument('--saveraw', action='store_true', help='Save sdssproc outputs')
    rungroup.add_argument('--debug', action='store_true', help='Save extraction debug files')
    rungroup.add_argument('--no_db', action='store_true', help='skip Database operations')
    rungroup.add_argument('--fast_no_db',  required=False,
                          help='When using --no_db, streamlines process and only gets parallax from MOS target files')
    rungroup.add_argument('--sdss_access_release', required=False,
                          help='sdss_access data release (defaults to sdsswork), required if you do not have proprietary access, '+
                               'otherwise see https://sdss-access.readthedocs.io/en/latest/auth.html#auth', default='sdsswork')
    rungroup.add_argument('--sdss_access_remote', help='allow for remote access to data using sdss-access', action='store_true')
    rungroup.add_argument('--dr19', help='Limit targeting flags to DR19 cartons', action='store_true')


    fieldgroup = parser.add_argument_group('Select Fields')
    fieldgroup.add_argument('--fieldids', '-f', nargs='*', help='Plate/Field numbers to reduce default="*"', type=str)
    fieldgroup.add_argument('--fieldstart', help='Starting Field/Plate number', default=None, type=str)
    fieldgroup.add_argument('--fieldend', help='End Field/Plate number', default=None, type=str)

    mjdgroup = parser.add_argument_group('Select MJDs')
    mjdgroup.add_argument('--mjd', '-m', nargs='*', help='MJD dates to reduce; default="*"', type=str)
    mjdgroup.add_argument('--mjdstart', help='Starting MJD', default=None, type=str)
    mjdgroup.add_argument('--mjdend', help='Ending MJD', default=None, type=str)

    slurmgroup = parser.add_argument_group('Slurm Options')
    slurmgroup.add_argument('--no_write', action='store_true', help='skip writing and submitting job')
    slurmgroup.add_argument('--kingspeak', action='store_true', help='Use kingspeak rather then notchpeak')
    slurmgroup.add_argument('--shared', action='store_true', help='Node sharing', default=share)
    slurmgroup.add_argument('--fast', action='store_true', help='Use SDSS fast queue')
    slurmgroup.add_argument('--mem_per_cpu', help='Memory allocated per CPU', type=str)
    slurmgroup.add_argument('--walltime', help='Wall time in hours', type=str)
    slurmgroup.add_argument('--nodes', default=1, help='Number of Nodes', type=int)
    slurmgroup.add_argument('--ppn', help='Number of processors per node', type=int)
    slurmgroup.add_argument('--nosubmit', action='store_true', help='Build, but not submit redux files')
    slurmgroup.add_argument('--clobber', action='store_true', help='Clobber redux')
    
    customgroup = parser.add_argument_group('Custom Coadd Options')
    customgroup.add_argument('--epoch', action='store_true', help = 'Epoch Coadds')
    customgroup.add_argument('--custom', help = 'Name of custom Coadd Schema', type=str)
    customgroup.add_argument('--allsky', action='store_true', help = 'All Sky Coadds')
    customgroup.add_argument('--coadd_only', action='store_true', help= 'Run spspec_target_merge only', dest = 'custom_coadd_only')
    customgroup.add_argument('--1dpost', action='store_true', help= 'Run 1d analysis and post processing only', dest='custom_1dpost')

    emailgroup = parser.add_argument_group('Email outputs')
    emailgroup.add_argument('--email', action='store_true', help='Email log using $HOME/daily/etc/emails')

    args = parser.parse_args()

    if args.bay15:
        args.map3d = 'bayestar15'
#    if args.eden23:
#        args.map3d = 'edenhofer2023'
    if args.merge3d:
        args.map3d = 'merge3d'
    if (args.sdssv is True) or (args.sdssv_fast is True) or (args.sdssv_noshare is True) :
        args.MWM_fluxer      = True
        args.no_reject       = True
        args.no_merge_spall  = True
        if args.mjdstart is None and args.mjd is None: 
            args.mjdstart    = 59146
        if args.fieldids is None and args.fieldstart is None:
            args.fieldstart  = 15000
        if args.walltime is None:
            args.walltime    = '40:00:00'
        if not args.sdssv_noshare:
            args.shared      = True
            if not args.kingspeak is True:
                if args.mem_per_cpu is None:
                    args.mem_per_cpu = '7500'
            else:
                if args.mem_per_cpu is None:
                    args.mem_per_cpu = '3750'
        else:
            if args.mem_per_cpu is None:
                args.mem_per_cpu = '125000'


    if args.sdssv_fast is True:
        args.fast            = True

    if args.walltime is None:
        args.walltime       = '40:00:00'
        if not args.kingspeak is True:
            if args.mem_per_cpu is None:
                args.mem_per_cpu = '7500'
        else:
            if args.mem_per_cpu is None:
                args.mem_per_cpu = '3750'

    if args.custom is not None:
        args.skip2d  = True

    args.include_bad = True
    args.xyfit       = True
    args.loaddesi    = True
                

    if args.lco is True:
        args.obs = ['lco']
    elif args.apo is True:
        args.obs = ['apo']
    if args.kingspeak is True:
        args.shared = False

    args_dic = vars(args)
    for key in ['sdssv','lco','apo','sdssv_fast', 'bay15', 'sdssv_noshare', 'merge3d']: #,'eden23']:
        args_dic.pop(key)

    queue = uubatchpbs(**args_dic)


