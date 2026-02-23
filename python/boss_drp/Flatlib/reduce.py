#!/usr/bin/env python3
from boss_drp.utils import (find_nearest_indx, load_env)
from boss_drp.field import Field
from boss_drp.prep.spplan_trace import spplanTrace

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
    warnings.warn('No slurm package installed: printing command to STDOUT for manual run',SlurmWarning)
    noslurm = True
    
from pydl.pydlutils.yanny import yanny, read_table_yanny
from astropy.io import fits
from astropy.table import Table
from os import getenv, symlink, makedirs
import os.path as ptt
from glob import glob
import argparse
import numpy as np
import sys
import io
from tqdm import tqdm
from shutil import copy2

def match_arc(arcs, flat, single_cal= True):
    m_arc = np.where(arcs['fieldid'].data == flat['fieldid'])[0]
    if len(m_arc) != 0:
        m_arc = m_arc[-1]
        arcs = Table(arcs[m_arc])
    else:
        if single_cal is True:
            texp = np.nanmean(flat['TAI'].data)
            try:
                idx = find_nearest_indx(arcs['TAI'].data, texp)
            except:
                return(None)
            arcs = Table(arcs[idx])
    return(arcs)
    
    
def create_run(dir_, specdir, mjd, obs='lco',no_run=False,
                submit=False, nodes =None, clobber=False, link=False):
    if type(mjd) is str: mjd=[mjd]
    
    cmds = []
    logs = []
    for mj in tqdm(mjd, desc='MJD',leave=False, position=0):
        legacy=False; plates=False
        if int(mj) < 59030: legacy= True
        elif int(mj) < 59550: plates= True
        idl = 'run_spcalib, mjd={mjd}'
        if plates is True: idl +=', /plates'
        if legacy is True: idl +=', /legacy'
        if obs.lower() == 'lco': idl +=', /lco'
        
        plan = read_table_yanny(ptt.join(dir_,'calibs',obs,str(mj),f'spPlanTrace-{mj}_{obs.upper()}.par'), 'SPEXP')
        plan.convert_bytestring_to_unicode()
        
        def get_tai(row, specdir, mjd):
            f = ptt.join(specdir,mjd,row['name'][0])
            try:
                hdr = fits.getheader(f)
            except:
                hdr = fits.getheader(f+'.gz')
            return(hdr['TAI-BEG'] + hdr['EXPTIME']/2.0)
            
        plan['TAI'] = [get_tai(row, specdir, mj) for row in plan]
        flats = np.where(plan['flavor'] == 'flat')[0]
        arcs = np.where(plan['flavor'] == 'arc')[0]
        outdir = ptt.join(dir_,'calibs',obs,str(mj))

        for f in tqdm(flats, desc='flat', leave=False, position=1):
            flatnames = plan['name'][f]
            field = plan['fieldid'][f]
            for i, flatname in enumerate(flatnames):
                flatinfoname = flatname.replace('.fit','.fits').replace('sdR','spFlat')
                if not clobber:
                    if ptt.exists(ptt.join(outdir, flatinfoname)):
                        continue
                if link:
                    fc = Field(getenv('BOSS_SPECTRO_REDUX'),getenv('RUN2D'),field)
                    ff = ptt.join(fc.dir(),flatinfoname)
                    if ptt.exists(ff):
                        if ptt.exists(ptt.join(outdir,flatinfoname)):
                            continue
                        symlink(ff,ptt.join(outdir,flatinfoname))
                        continue
                    elif ptt.exists(ff+'.gz'):
                        if ptt.exists(ptt.join(outdir,flatinfoname+'.gz')):
                            continue
                        symlink(ff+'.gz',ptt.join(outdir,flatinfoname)+'.gz')
                        continue
                if no_run:
                    continue
                plottitle = flatinfoname.replace('.fits','.ps')
                marcs = match_arc(plan[arcs], plan[f])
                if marcs is None:
                    continue
                for arc in marcs:
                    arcname = arc['name'][i]
                    arcinfoname = arcname.replace('.fit','fits').replace('sdr','spArc')
                    try:
                        cartid = fits.getheader(ptt.join(specdir,f'{mj}',flatname))['cartid']
                    except:
                        cartid = fits.getheader(ptt.join(specdir,f'{mj}',flatname+'.gz'))['cartid']
                    idl = ('spcalib, {flatname}, {arcname}, cartid={cartid}, rawdir={specdir}, '+
                                     'plottitle={plottitle}, timesep=0, outdir={outdir}, '+
                                     'flatinfoname={flatinfoname}, arcinfoname={arcinfoname} ')
                    if int(mj)< 59030:
                        idl += ', /legacy'
                    elif int(mj)<59550:
                        idl += ', /plates'
                    idl = idl.format(flatname=flatname,arcname=arcname,cartid=cartid,
                                     specdir=specdir,plottitle=plottitle,outdir=outdir,
                                     flatinfoname=flatinfoname, arcinfoname=arcinfoname)
                    script = 'idl -e "{idl}"'.format(idl=idl)
                    logfile = ptt.join(dir_,'logs',f'{mjd}', flatinfoname.replace('.fits','.log'))
                    cmds.append(script)
                    logs.append(logfile)
            
    if no_run:
        return
    if len(cmds)> 0:
        return
    
    if not noslurm:
        queue1 = queue(verbose=True)
        

    share=True
    maxnodes = 5
    alloc = getenv('SLURM_ALLOC')
    maxcore = int(getenv('SLURM_PPN'))
    if nodes is not None:
        maxnodes = Nones
    ncmds = len(cmds)
    nnodes = int(np.ceil(ncmds/maxcore))
    
    if nnodes > maxnodes:
        nnodes=maxnodes
    ncore = int(np.ceil(nmjds/nnodes))
    if ncore > maxcore:
        ncore = maxcore
    
        
    if len(mjd) == 1: mjs = ' '+mjd[0]
    else: mjs=''

    if not noslurm:
        queue1.create(label='BOSSFlatlib_'+getenv('RUN2D')+mjs, nodes=nnodes,
                    ppn=ncore, qos='sdss',  shared=share, walltime='168:00:00', alloc=alloc)

    for i, cmd in enumerate(cmds):
        if not noslurm:
            queue1.append(cmd, outfile = logs[i]+'.o', errfile = logs[i]+'.e')
        else:
            print(cmd+' > '+logs[i]+'.o 2> '+logs[i]+'.e')
    queue1.commit(hard=True, submit=submit)


def reduce(dir_, mjd, link=False, lco=False, plates=False, nodes=None,no_run=False,
                  legacy=False, fps=False, nosubmit=False, deep=False, link_all=False,
                  link_traceflat=False, mjdstart = None):

    if fps or legacy or plates:
        specdir= 'BOSS_SPECTRO_DATA_S' if lco else 'BOSS_SPECTRO_DATA_N'
        specdir = getenv(specdir)
        mjd=[ptt.basename(x) for x in glob(ptt.join(specdir,'?????'))]
        mjd=np.asarray(mjd,dtype=int)
        mjd = mjd[np.where(mjd >= 59550)[0]] if fps else mjd[np.where(mjd < 59550)[0]]
        if mjdstart is not None:
            mjd = mjd[np.where(mjd >= mjdstart)[0]]
        if plates:
            mjd = mjd[np.where(mjd >= 59030)[0]]
        if legacy:
            mjd = mjd[np.where(mjd <59030)[0]]
        mjd = mjd.astype(str).tolist()
        
    obs = 'lco' if lco else 'apo'
    makedirs(ptt.join(dir_,'calibs',obs), exist_ok = True)
    if not deep:
        libmjd = [ptt.basename(ptt.dirname(x)) for x in glob(ptt.join(dir_,'calibs',obs,'?????','spPlanTrace*'))]
    else:
        libmjd = []
    mjds=[]
    for mj in mjd:
        if mj not in libmjd:
            if ptt.exists(ptt.join(getenv('BOSS_SPECTRO_REDUX'),getenv('RUN2D'),
                         'trace',mj,f'spPlanTrace-{mj}_{obs.upper()}.par')):
                makedirs(ptt.join(dir_,'calibs',obs,mj), exist_ok= True)
                if not ptt.exists(ptt.join(dir_,'calibs', obs,mj,
                            f'spPlanTrace-{mj}_{obs.upper()}.par')):
                    copy2(ptt.join(getenv('BOSS_SPECTRO_REDUX'),getenv('RUN2D'),
                                   'trace',mj,f'spPlanTrace-{mj}_{obs.upper()}.par'),
                          ptt.join(dir_,'calibs', obs,mj,
                                f'spPlanTrace-{mj}_{obs.upper()}.par'))
                mjds.append(mj)
                continue
            if not no_run:
                makedirs(ptt.join(dir_,'calibs',obs,mj), exist_ok= True)

                nmjds = spplanTrace(mjd=mj,lco=lco, legacy=legacy, plates=plates,
                                   sav_dir=ptt.join(dir_,'calibs',obs), flib=True)
            
                if nmjds is not None:
                    mjds.append(mj)

    if link_all:
        sp = '2' if lco else '1'
        for f in tqdm(glob(ptt.join(getenv('BOSS_SPECTRO_REDUX'),getenv('RUN2D'),
                            f'??????','??????',f'spFlat-?{sp}-????????.fits.gz'))):
            tmjd = str(fits.getval(f,'MJD'))
            if tmjd in mjd:
                makedirs(ptt.join(dir_,'calibs',obs,tmjd), exist_ok=True)
                if not ptt.exists(ptt.join(dir_,'calibs',obs,tmjd,ptt.basename(f))):
                    symlink(ptt.abspath(f), ptt.join(dir_,'calibs',obs,tmjd,ptt.basename(f)))
    if link_traceflat:
        sp = '2' if lco else '1'
        for f in tqdm(glob(ptt.join(getenv('BOSS_SPECTRO_REDUX'),getenv('RUN2D'),
                            'trace','?????',f'spTraceFlat-?{sp}-????????.fits.gz'))):
            tmjd = str(fits.getval(f,'MJD'))
            if tmjd in mjd:
                makedirs(ptt.join(dir_,'calibs',obs,tmjd), exist_ok=True)
                if not ptt.exists(ptt.join(dir_,'calibs',obs,tmjd,ptt.basename(f))):
                    symlink(ptt.abspath(f), ptt.join(dir_,'calibs',obs,tmjd,ptt.basename(f.replace('spTraceFlat','spFlat'))))
        
    create_run(dir_, specdir, mjds, obs = obs, submit=(not nosubmit),
                link=link, nodes =nodes, no_run=no_run)

