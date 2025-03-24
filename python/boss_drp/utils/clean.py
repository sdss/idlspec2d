#!/usr/bin/env python3
from boss_drp.field import *

from collections import OrderedDict
from pydl.pydlutils import yanny
import os.path as ptt
from os import getcwd, getenv, rmdir
import numpy as np
from glob import glob
from os import remove as rm
import sys
import subprocess
from astropy.table import Table, vstack
from tqdm import tqdm

def remove(pt, dry= True):
    for file in glob(pt):
        try:
            if ptt.exists(file):
                if not dry:
                    if ptt.isdir(file):
                        rmdir(file)
                    else:
                        rm(file)
                else:
                    tqdm.write(f'\t {dry} {file}')
        except Exception as e:
            tqdm.write(f'{e}')
            pass

fmag_f = ['spfibermap-{field}-{mjd}.fits','spfibermap-{field}-{mjd}.log']

debug_f = ['*.inp','*_supp.fits','spProc_bksub-*-{expid}.fits','sdProc-*{expid}.fits',
           'extraction/extract_*{expid}*.prt','extraction/gauss_*{expid}*.prt',
           'extraction/sigma_*{expid}*.prt',
           'flat_extraction/extract_*{expid}*.prt','flat_extraction/gauss_*{expid}*.prt',
           'spFrame_bksub-*{expid}+.fits','fitspectraresol-*{expid}*',
           'spArcFlux-*{expid}.fits','spProc-*{expid}.fits','lamps.html',
           'spFrame_bksub-*{expid}.fits','spFrame_preext-*{expid}.fits',
           'spFrame_preflat-*{expid}.fits','extraction','flat_extraction']
spec2d_f = ['spec2d-{field}-{mjd}.started','spec2d-{field}-{mjd}.done',
            'spDiag2d-{field}-{mjd}.log*','spDiag2d-{field}-{mjd}.ps*',
            'spDiag2d-{field}-{mjd}.pdf*',
            'spArc-*{expid}.fits.gz','spFluxdistort-*{expid}.fits',
            'spFlat_*{expid}.fits.gz','spFlat-*{expid}.fits.gz',
            'spFlat_*{expid}.fits', 'spFlatFlux-*{expid}.fits',
            'spFrame-*-{expid}.fits.gz']

comb_f = ['specombine-{field}-{mjd}.started','specombine-{field}-{mjd}.done']
comb_f_epoch = ['spPlancombepoch-{field}-{mjd}.started','spPlancombepoch-{field}-{mjd}.done']
comb_f_ext =  ['spDiagcomb-{field}-{mjd}.log*','spDiagcomb-{field}-{mjd}.ps*',
               'spDiagcomb-{field}-{mjd}.pdf*',
               'spCFrame-*{expid}.fits.gz','spCFrame-*{expid}.fits',
               'spFluxcorr-*{expid}.fits*','spFluxcalib-*{expid}.fits*',
               'spField-{field}-{mjd}.fits','coadd/{mjd}/spSpec-{field}-{mjd}-*.fits',
               'coadd/{mjd}','coadd','spSN2d-{field}-{mjd}-*.ps*',
               'spSN2d-{field}-{mjd}-*.pdf*','spSN2d-{field}-{mjd}.pdf*',
               'spSN2d-{field}-{mjd}.ps*','spFluxdistort-*{expid}.fits*',
               'spFluxdistort-*{expid}.ps*','spFluxdistort-*{expid}.pdf*',
               'spFluxdistort-{field}-{mjd}.pdf*','spFluxdistort-{field}-{mjd}.ps*']

spec1d_f = ['spec1d-{field}-{mjd}.started','spec1d-{field}-{mjd}.done',
            '{run1d}/spDiag1d-{field}-{mjd}.log*',
            '{run1d}/spDiag1d-{field}-{mjd}.ps*',
            '{run1d}/spDiag1d-{field}-{mjd}.pdf*',
            '{run1d}/spZall-{field}-{mjd}.fits',
            '{run1d}/spZbest-{field}-{mjd}.fits',
            '{run1d}/spZline-{field}-{mjd}.fits',
            '{run1d}/spXCSAO-{field}-{mjd}.fits',
            '{run1d}/spXCSAO-{field}-{mjd}.log*','{run1d}']

post_f = ['fieldlist-{field}-{mjd}.log*']
merge_f = ['fieldmerge-{field}.log*','fieldmerge-{field}-{mjd}.log*',
          'spAll-{field}-{mjd}.log*']
specF_log_f = ['spSpec_reformat-{field}-{mjd}.log*','spec-{field}-{mjd}.log*']
specF_f = ['spec-{field}-{mjd}-*.fits','spAll-{field}-{mjd}.fits.gz',
           'spAllLine-{field}-{mjd}.fits.gz']
specI_f = ['*','']
calib_f = ['spCalib_QA-{run2d}-{field}-{mjd}.ps*','spCalib_QA-{run2d}-{field}-{mjd}.log*',
           'spCalib_QA-{run2d}-{field}-{mjd}.pdf*']
reset_f = ['spPlan2d-{field}-{mjd}.par','spPlancomb-{field}-{mjd}.par']
reset_f_epoch = ['spPlancombepoch-{field}-{mjd}.par']
redux_f = ['redux-{field}-{mjd}*']

def clean_fmjd(topdir, run2d, run1d, field, mjd, epoch=False, dry=False,
               clean_type='full', reset=False, remove_redux = False,
               verbose=True):
    field = field_to_string(field)
    if not epoch:
        es = ''
        plan =  f'spPlancomb-{field}-{mjd}.par'
    else:
        es = 'epoch'
        plan = f'spPlancombepoch-{field}-{mjd}.par'
    fc = Field(topdir,run2d, field)
    plan1d = yanny.read_table_yanny(ptt.join(fc.dir(),es, plan),'SPEXP')
    spPlan2ds = plan1d.meta['planfile2d'].replace("'",'').split()
    sp2d = Table()
    for spPlan2d in spPlan2ds:
        plan2d = yanny.read_table_yanny(ptt.join(fc.dir(),spPlan2d),'SPEXP')
        plan2d.meta = {}
        sp2d = vstack([sp2d, plan2d])
    sp2d['expid'] = ''
    sp2d['expid'] = [ptt.splitext(x)[0].split('-')[-1] for x in np.asarray(sp2d['name'][:,0],dtype=str).tolist()]
    sp1d = plan1d
    sp1d['expid'] = ''
    sp1d['expid'] = [ptt.splitext(x)[0].split('-')[-1] for x in np.asarray(sp1d['name'][:,0],dtype=str).tolist()]
        
    if not epoch:
        tcomb_f = comb_f + comb_f_ext
    else:
        tcomb_f = comb_f_epoch +comb_f_ext
    
    if clean_type == 'all':
        clean = OrderedDict({'fmag':fmag_f, 'debug':debug_f,'spec2d':spec2d_f,'comb':tcomb_f,
                 'spec1d':spec1d_f,'post':post_f,'merge':merge_f,'specF_log':specF_log_f,
                 'specF':specF_f,'specI':specI_f,'spCalib':calib_f})
        if reset:
            remove_redux = True
            if epoch:
                clean['reset'] = reset_f_epoch
            else:
                clean['reset'] = reset_f
        if remove_redux:
            clean['redux'] = redux_f
    elif clean_type == 'spec2d':
        clean = OrderedDict({'debug':debug_f,'spec2d':spec2d_f,
                 'comb':tcomb_f,'spec1d':spec1d_f,
                 'post':post_f,'merge':merge_f,'specF_log':specF_log_f,
                 'specF':specF_f,'specI':specI_f,'spCalib':calib_f})
    elif clean_type == 'comb':
        clean = OrderedDict({'comb':tcomb_f,'spec1d':spec1d_f,
                 'post':post_f,'merge':merge_f,'specF_log':specF_log_f,
                 'specF':specF_f,'specI':specI_f,'spCalib':calib_f})
    elif clean_type == 'spec1d':
        clean = OrderedDict({'spec1d':spec1d_f, 'post':post_f,'merge':merge_f,
                 'specF_log':specF_log_f, 'specF':specF_f,'specI':specI_f,'spCalib':calib_f})
    elif clean_type == 'post':
        clean = OrderedDict({'post':post_f,'merge':merge_f,'specF_log':specF_log_f,
                 'specF':specF_f,'specI':specI_f,'spCalib':calib_f})
    elif clean_type == 'merge':
        celan = OrderedDict({'merge':merge_f,'specF_log':specF_log_f,
                 'specF':specF_f,'specI':specI_f,'spCalib':calib_f})
    elif clean_type == 'reformat':
        clean = OrderedDict({'specF_log':specF_log_f,'specF':specF_f,
                'specI':specI_f, 'spCalib':calib_f})
    elif clean_type == 'spCalib':
        clean = OrderedDict({'spCalib':calib_f})
    else:
        exit()


    for key, value in clean.items():
        if key in ['spec2d', 'debug']:
            expTab = sp2d
        elif key in ['comb']:
            expTab = sp1d
        else:
            expTab = None
        _step_clean(key,value, topdir, run2d, run2d, field, mjd, expTab,
                    epoch=epoch, dry=dry, verbose=verbose)

def _step_clean(step, fpath, topdir, run2d, run1d, field, mjd, expTab,
                epoch=False, dry=False, verbose = True):
    es = '' if not epoch else 'epoch'
    fc = Field(topdir, run2d, field, epoch = epoch)
    if step not in ['specI','specF']:
        basedir = fc.dir()
    elif step == 'specI':
        basedir = fc.png_dir(run1d, mjd)
    elif step == 'specF':
        basedir = fc.spec_dir(mjd)
    for fp in fpath:
        _type_clean(basedir, fp, {'run2d':run2d,'run1d':run1d,'field':field,'mjd':mjd},
                    expTab, dry=dry, verbose=verbose)
        if step == 'specF':
            _type_clean(basedir.replace('full','lite'), fp,
                        {'run2d':run2d,'run1d':run1d,'field':field,'mjd':mjd}, expTab,
                        dry=dry, verbose=verbose)

def _type_clean(basedir, fpath, params, expTab, dry=False, verbose=True):
    if 'expid' not in fpath:
        expids = ['']
    else:
        expids = []
        for row in expTab:
            expids.append(row['expid'])
    for expid in expids:
        if verbose:
            tqdm.write(ptt.join(basedir,fpath.format(**params,expid=expid)))
        remove(ptt.join(basedir,fpath.format(**params,expid=expid)), dry=dry)
    

