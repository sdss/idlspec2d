#!/usr/bin/env python3
from boss_drp.utils import get_dirs, load_env, mjd_match
from boss_drp.field import Field, field_to_string
from boss_drp.utils import jdate
from boss_drp import daily_dir
from boss_drp.utils.splog import splog

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
import argparse
from os import getenv
import os.path as ptt
from glob import glob
from pydl.pydlutils.yanny import yanny, read_table_yanny
import sys
import numpy as np
import io
import datetime
import astropy.time
import pandas as pd

if getenv('SLURM_VERS') == 'notchpeak': 
    share = True
else: 
    share = False


def build_cmd(topdir=None,run2d=None,run1d=None,idlutils_1d=None,
        no_reject=False, MWM_fluxer=False, map3d='bayestar15',
        noxcsao=False, skip_specprimary=False,
        no_merge_spall=False, skip2d=False, onestep_coadd=False,
        fibermap_clobber=False, lco=False, plan2d=None, plancomb=None,
        fieldmjd=None, post_idl=False, only1d=False, daily=False,
        custom=None, allsky=False, epoch=False, saveraw=False, debug=False,
        release = None, remote = False, force_arc2flat=False, nodist=False,
        no_db=False, fast_no_db=False,no_healpix=False, dr19=False,
        custom_coadd_only=False, custom_1dpost=False, redux=None, a2t=False, **kwargs):

    field = fieldmjd.split('-')[-2]
    mjd = fieldmjd.split('-')[-1]

    cmd = []
    cmd.append('# Auto-generated batch file '+datetime.datetime.now().strftime("%c"))

    fc = Field(topdir, run2d, field, epoch=epoch, custom_name=custom, custom=allsky)
    cmd.append('cd '+fc.dir())
    if custom is not None:
        fieldmjd = field+'-'+mjd

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
    spSpecRef_key = "" # " --lsdr10"
    flist_key = ""
    fmerge_key = " --include_bad"
    if run2d.lower() == 'master':
        spSpecRef_key =" --lsdr10"

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
    if release is not None:
        fmap_keywords = fmap_keywords + " --release "+release
    if remote is not None:
        fmap_keywords = fmap_keywords + " --remote"

    if lco:
        spreduce2d_keys = spreduce2d_keys + ' /lco,'
    if MWM_fluxer:
        spreduce2d_keys = spreduce2d_keys + ' /MWM_fluxer,'
        rm_combine_keys = rm_combine_keys + ' /MWM_fluxer,'
        if map3d:
            spreduce2d_keys = spreduce2d_keys + f' map3d={map3d},'
    if saveraw:
        spreduce2d_keys = spreduce2d_keys + ' /saveraw,'
    if debug:
        spreduce2d_keys = spreduce2d_keys + ' /debug,'
    
    lco_rm_fps = [23129,23130,23131,23132,23133,23134,23135,23136,23137,23175,
                  23288,23408,23409,23410,31687,31688,112357,112358,112362]
    if lco and int(mjd) < 60402 and int(field) not in lco_rm_fps:
            spreduce2d_keys = spreduce2d_keys + ' /force_arc2trace,'
    if no_reject:
        rm_combine_keys = rm_combine_keys + ' /no_reject,'
    if onestep_coadd:
        rm_combine_keys = rm_combine_keys + ' /onestep_coadd,'
    if nodist:
        rm_combine_keys = rm_combine_keys + ' /nodist,'
    if a2t:
        spreduce2d_keys = spreduce2d_keys + ' /force_arc2trace,'
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
        fmerge_key = fmerge_key + pycustomkey
        flist_key = flist_key + pycustomkey
        spSpecRef_key = spSpecRef_key + pycustomkey
        xcsao_keys = xcsao_keys       + pycustomkey
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
                cmd.append(f"readfibermaps --spplan2d {plan}")
            else:
                cmd.append(f"readfibermaps --spplan2d {plan} --dr19")
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
            cmd.append("echo 'rm_combine_script, "+'"spPlancombepoch-'+fieldmjd+'.par",'+
                        rm_combine_keys+' run2d="'+run2d+'"'+"' |idl")
            cmd.append('touch spPlancombepoch-'+fieldmjd+'.done')
        elif custom is None:
            cmd.append('touch specombine-'+fieldmjd+'.started')
            cmd.append("echo 'rm_combine_script, "+'"spPlancomb-'+fieldmjd+'.par",'+
                        rm_combine_keys+' run2d="'+run2d+'"'+"' |idl")
            cmd.append('touch specombine-'+fieldmjd+'.done')


    if run1d is None: 
        run1d=run2d
    if run2d != run1d: 
        cmd.append('module switch idlspec2d idlspec2d/'+run1d)
    if idlutils_1d is not None:
        cmd.append('module switch idlutils idlutils/'+idlutils_1d)

    if allsky is False:
        cmd.append('touch spec1d-'+fieldmjd+'.started')
        cmd.append("echo 'spreduce1d_empca, "+'"spField-'+fieldmjd+'.fits",'+
                    spec1d_keys+' run1d="'+run1d+'"'+"' |idl")
        if not noxcsao:
            cmd.append('run_PyXCSAO spField-'+fieldmjd+'.fits'+xcsao_keys+' --run1d "'+run1d+'"')
        cmd.append('touch spec1d-'+fieldmjd+'.done')
    else:
        global plancomb_last
        global EPOCH_COMBINE
        global EPOCH_OBS
        if 'plancomb_last' not in globals():
            plancomb_last = ''
        if plancomb != plancomb_last:
            cplan = read_table_yanny(plancomb, 'COADDPLAN')
            EPOCH_COMBINE = np.sort(np.unique(cplan['EPOCH_COMBINE']))
            try:
                EPOCH_OBS = cplan.meta['OBS'].lower()
            except:
                EPOCH_OBS = None
            plancomb_last = plancomb
        for mjec in EPOCH_COMBINE:
            if custom_1dpost:
                test_mjec = redux.split('_')[-1]
                if str(mjec) != test_mjec:
                    continue
            if EPOCH_OBS is not None:
                fmjd = custom+'_'+EPOCH_OBS+'-'+str(mjec)
            else:
                fmjd = custom+'-'+str(mjec)
            cmd.append('touch spec1d-'+fmjd+'.started')
            cmd.append("echo 'spreduce1d_empca, "+'"spFullsky-'+fmjd+'.fits",'+
                        spec1d_keys+' run1d="'+run1d+'"'+"' |idl")
            if not noxcsao:
                cmd.append('run_PyXCSAO spFullsky-'+fmjd+'.fits'+xcsao_keys+' --run1d "'+run1d+'"')
            cmd.append('touch spec1d-'+fmjd+'.done')


    if not only1d:
        if allsky is False:
            fmjd = [fieldmjd]
        else:
            fmjd = []
            for mjec in EPOCH_COMBINE:
                if custom_1dpost:
                    test_mjec = redux.split('_')[-1]
                    if str(mjec) != test_mjec:
                        continue


                if EPOCH_OBS is not None:
                    fmjd.append(custom+'_'+EPOCH_OBS+'-'+str(mjec))
                else:
                    fmjd.append(custom+'-'+str(mjec))
                
        for fm in fmjd:
            field, mjd = fm.split('-')
            if allsky is False:
                cmd.append('')
                cmd.append('#- Make Field List')
                cmd.append(f"fieldlist --create --run1d {run1d} --run2d {run2d} {flist_key} --logfile fieldlist-{field}-{mjd}.log")
            
                cmd.append('')
                cmd.append('#- Make Field spAll file')
                cmd.append(f"fieldmerge {fmerge_key} --clobber --run2d {run2d} --field {field} --mjd {mjd}")
                cmd.append('')
                cmd.append('#- Make final spectra files')
                cmd.append(f"spSpec_reformat --field {field} --mjd {mjd} -p --run1d {run1d} --run2d {run2d} {spSpecRef_key}")
            else:
                cmd.append('')
                cmd.append('#- Make Field spAll file')
                cmd.append(f"fieldmerge {fmerge_key} --clobber --run2d {run2d} --field {field} --custom {custom} --run1d {run1d} --mjd {mjd}")
                cmd.append('')
                cmd.append('#- Make final spectra files')
                cmd.append(f"spSpec_reformat --field {field} --mjd {mjd} -p --run1d {run1d} --run2d {run2d} --custom {custom} {spSpecRef_key} ")

            if allsky is False:
                cmd.append('')
                cmd.append('#- Make Spectro-Photometric Calibration QA Figures')
                cmd.append("echo 'spcalib_qa, fieldid="+field+", mjd="+mjd+','+spcalib_keywords+' run2d="'+run2d+'"'+"' |idl")

        if not no_merge_spall:
            cmd.append('')
            cmd.append('#- Making Summary Files')
            if allsky is False:
                cmd.append(f"fieldlist --create --run1d {run1d} --run2d {run2d} {flist_key}")
                cmd.append(f"fieldmerge --lite {fmerge_key} --run2d {run2d} --remerge_fmjd {field}-{mjd} --update_specprimary ")
            else:
                cmd.append(f"fieldmerge --lite {fmerge_key} --custom {field} --run1d {run1d} --run2d {run2d} --remerge_fmjd {field}-{mjd} --update_specprimary")

                    
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
               fibermap_clobber = False, nbundle = None,
               only1d = False, force_arc2flat=False,
               field = None, fieldstart = None, fieldend = None,
               mjd = None, mjdstart = None, mjdend = None, saveraw=False,
               no_write = False, shared = share, fast = False, mem = None, 
               mem_per_cpu = getenv('SLURM_MEM_PER_CPU'), walltime = '336:00:00',
               nodes = None, ppn = None, nosubmit = False, daily=False,
               debug=False, no_db=False, dr19=False, nodist=False, 
               clobber = False, custom= None, allsky = False, epoch = False, no_healpix=False,
               email = False, logger=None, fast_no_db=False,
               remote = False, release=None,allemail=False,
               custom_coadd_only=False, custom_1dpost=False, a2t=False):


    if daily is False:
        splog.emailer()
    
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
    if cmdinputs['custom'] is None:
        cmdinputs.pop('allsky')
    for key in list(cmdinputs.keys()):
        if cmdinputs[key] is None:
            cmdinputs.pop(key)
    error = False
    
    if not ptt.isdir(topdir):
        splog.warning('topdir (BOSS_SPECTRO_REDUX) is invalid')
        error = True
    else:
        splog.info('topdir   '+topdir)
    splog.info(pd.Series(cmdinputs).to_string())
    
    afc = Field(topdir, run2d, '*')
    if not allsky:
        fielddirs = get_dirs(ptt.dirname(afc.dir()), field=True,
                             match = field, start=fieldstart, end=fieldend)
        if len(fielddirs) == 0:
            splog.warning('No Directories Found')
            error = True
    else:
        if len(obs) == 2:
            fielddirs = [custom]
        else:
            fielddirs = []
        fielddirs.extend([custom+'_'+ob for ob in obs])
        
    if error:
        splog.close_elogger()
        return()

    redux_list = []
    skipped = 0
   
    if epoch is False:
        if custom is None:
            plan_str = 'spPlancomb-*.par'
        else:
            plan_str = 'spPlanCustom-{custom}-*.par'
            plan_str_bkup = plan_str
    else:
        if custom is None:
            plan_str = 'spPlancombepoch-*.par'
        else:
            plan_str = 'spPlancombepoch_{custom}-*.par'
            plan_str_bkup = plan_str

    if clobber:
        if mjd is None and mjdstart is None and mjdend is None:
            splog.info('No MJDs Selected while clobber is set')
            val = input('Do you want to continue? (yes/NO)')
            if val.lower() != 'yes':
                exit()

    for fielddir in fielddirs:
        if custom is not None:
            plan_str = plan_str_bkup
            plan_str = plan_str.format(custom=fielddir)
        fc = Field(topdir, run2d, fielddir, custom_name = custom, epoch = epoch)
        planfile = glob(ptt.join(fc.dir(), plan_str))
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
            if thisobs.lower() == 'lco':
                lco = True
            else:
                lco = False

            if custom is None:
                fieldmjd = hdr['fieldid']+'-'+hdr['MJD']
                redux = ptt.join(fc.dir(), 'redux-'+fieldmjd)
            else:
                if allsky:
                    fieldmjd = fielddir+'-'+hdr['CreateMJD']
                else:
                    fieldmjd = fielddir+'-'+hdr['fieldid']+'-'+hdr['CreateMJD']
                if not custom_1dpost:
                    redux = ptt.join(fc.dir(), 'redux_'+fieldmjd)
                else:
                    redux = []
                    cplan = read_table_yanny(plan, 'COADDPLAN')
                    EPOCH_COMBINE = np.sort(np.unique(cplan['EPOCH_COMBINE']))
                    for mjec in EPOCH_COMBINE:
                        redux.append(ptt.abspath(ptt.join(fc.dir(),'redux_'+fieldmjd+'_'+str(mjec))))
                        if epoch: print(redux[-1])

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
    splog.info('')
    splog.info('---------------------------------------------------')
    splog.info('boss_redux: #BOSS Field-MJDs Done  = '+str(skipped))
    splog.info('            #BOSS Field-MJDs To Do = '+str(len(redux_list)))
    splog.info('---------------------------------------------------')
    splog.info('')
    
    if len(redux_list) == 0: 
        splog.close_elogger()
        return None, redux_list


    qos = 'sdss'
    alloc = getenv('SLURM_ALLOC')
    if fast is True: alloc = alloc+'-fast'
    partition = 'sdss-np'
    max_c = int(getenv('SLURM_PPN'))
    if ppn is None: ppn= np.min([max_c, np.max([len(redux_list),2])])


    if daily is True:
        if nodes is None and len(redux_list) > ppn:
            nodes = int(np.ceil(len(redux_list)/max_c))
        else:
            nodes = 1
    elif nodes is None:
        nodes = 1
    
        
    obsstr = '_'.join(obs).upper()
    mjdstr = str(mjd[0]) if daily else 'batch'
    if custom is None:
        if not epoch:
            label = f'{run2d}/{obsstr}/{mjdstr}/daily/'
        else:
            label = f'{run2d}/{obsstr}/{mjdstr}/epoch/'
    else:
        if not epoch:
            label = f'{run2d}/{obsstr}/{mjdstr}/{custom}/'
        else:
            label = f'{run2d}/{obsstr}/{mjdstr}/epoch/{custom}/'


    if nodes > 1:
        label = label.replace('/','_')

    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    if not no_write:
        queue1 = queue(key=None, verbose=True)
        bundle = True if nbundle is not None else False
        queue1.create(label=label, nodes=str(nodes), ppn=str(ppn),
                      partition = partition, alloc=alloc, shared=shared,
                      walltime=walltime, mem_per_cpu=mem_per_cpu, mem=mem,
                      nbundle = nbundle, bundle = bundle,)
    else:
        queue1 = None
    rlist = redux_list.copy()
    for i in range(len(redux_list)):
        cmd, log, err = make_run_cmd(redux_list[0])
        redux_list.pop(0)
        if not no_write:
            queue1.append(cmd, outfile = log, errfile = err)
        else: 
            splog.info(f'{cmd}  > {log} 2> {err}')
    if not no_write:
        queue1.commit(hard=True, submit=(not nosubmit))
    output = new_stdout.getvalue()
    sys.stdout = old_stdout
    splog.info(output)

    if email is True:
        if daily:
            tjdate = mjd[0]
        else:
            tjdate = jdate.astype(str)
        splog.send_email('UUBATCH '+run2d +' MJD='+str(tjdate) +' OBS='+','.join(obs),
                              ptt.join(daily_dir, 'etc','emails'), allemail=allemail)
    else:
        splog.close_elogger()

    return(queue1, rlist)

def make_run_cmd(redux):
    cmd = 'source '+redux
    log = redux+'.o'
    err = redux+'.e'
    return(cmd, log, err)

