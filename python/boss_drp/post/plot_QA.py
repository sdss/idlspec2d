#!/usr/bin/env python3
from boss_drp.field import field_spec_dir
from boss_drp.utils import load_env
from boss_drp import daily_dir, favicon, idlspec2d_dir
from boss_drp.utils import match as wwhere


try:
    from sdssdb.peewee.sdss5db.targetdb import database
    test = database.set_profile(load_env('DATABASE_PROFILE', default='pipelines'))
    from sdssdb.peewee.sdss5db.targetdb import Field, Cadence, DesignMode, Design, DesignToField
    from sdssdb.peewee.sdss5db import opsdb
    from sdssdb.peewee.sdss5db.opsdb import Exposure, CameraFrame, Camera
except:
    if load_env('DATABASE_PROFILE', default='pipelines').lower() in ['pipelines','operations']:
        print('ERROR: No SDSSDB access')
        exit()
    else:
        print('WARNING: No SDSSDB access')

from pydl.pydlutils.yanny import read_table_yanny

from astropy.io import fits
from astropy.table import Table, unique, vstack
from astropy.time import Time

import pandas as pd
from os import getenv, environ, makedirs
import numpy as np
import argparse
import os.path as ptt
import time
import logging
import warnings
import datetime
import re
from tqdm import tqdm
from glob import glob



import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.colors as pc
    colors={'g':'green','y':'goldenrod','r':'red','m':'magenta','b':'blue', 'k':'black'}
    axopts = dict( gridcolor='lightgrey', linecolor='darkgrey', showline=True,
                   mirror=True, tickcolor='darkgrey', ticks='outside',tickformat='.0f',
                   tickprefix='',ticksuffix='',tickmode='auto',
                   tickformatstops=[dict(dtickrange=[0,10], value='.2f'),
                                    dict(dtickrange=[10,None],value='.0f')])
    ayopts = dict( gridcolor='lightgrey', linecolor='darkgrey', showline=True,
                   mirror=True, tickcolor='darkgrey', ticks='outside',tickmode='auto')
    import plotly
    from jinja2 import Template
except:
    plotly = None

filters = ['G','R','I']
class Formatter(logging.Formatter):
    def __init__(self):
        super().__init__(fmt="%(levelno)d: %(msg)s", datefmt=None, style='%')
    def format(self, record):
        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        elif record.levelno == logging.DEBUG:
            self._style._fmt = '%(funcName)s: %(message)s'
        else:
            self._style._fmt = "%(levelname)s: %(message)s"
        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)
        # Restore the original format configured by the user
        self._style._fmt = format_orig
        return result


def load_fields(clobber_lists=False):
    if not ptt.exists(ptt.join(daily_dir,'etc', 'RM_fields')) or clobber_lists is True:
        try:
            rm_fields = []
            field = DesignToField.select().join(Design).join(DesignMode).where(DesignMode.label == 'dark_rm').switch(DesignToField).join(Field)
            for t in field: rm_fields.append(t.field.field_id)
            field = DesignToField.select().join(Design).join(DesignMode).where(DesignMode.label == 'dark_rm_eng').switch(DesignToField).join(Field)
            for t in field: rm_fields.append(t.field.field_id)
            rm_fields=list(set(rm_fields))
            trm_fields=np.char.add(np.asarray(rm_fields).astype(str),'\n').tolist()

            if ptt.exists(ptt.join(daily_dir,'etc')):
                with open(ptt.join(daily_dir,'etc', 'RM_fields'), 'w') as f:
                    f.writelines(trm_fields)
        except:
            rm_fields = []
            if not cron: warnings.warn('RM_fields file is missing and DB query failed', UserWarning)
            else: splog.warnings('Warning: RM_fields file is missing and DB query failed')
    else:
        with open(ptt.join(daily_dir,'etc', 'RM_fields'), 'r') as f:
            rm_fields = np.array(f.readlines())
            rm_fields= np.char.replace(rm_fields, '\n','').astype(int).tolist()

    if ptt.exists(ptt.join(daily_dir,'etc', 'RM_plates')):
        with open(ptt.join(daily_dir,'etc', 'RM_plates')) as f:
            rm_plates = np.array(f.readlines())
            rm_plates= np.char.replace(rm_plates, '\n','').astype(int).tolist()
            rm_fields.extend(rm_plates)

    if not ptt.exists(ptt.join(daily_dir,'etc', 'DarkMonitor_fields')) or clobber_lists is True:
        try:
            monit_fields = []
            field = DesignToField.select().join(Design).join(DesignMode).where(DesignMode.label == 'dark_monit').switch(DesignToField).join(Field)
            for t in field: monit_fields.append(t.field.field_id)
            field = DesignToField.select().join(Design).join(DesignMode).where(DesignMode.label == 'dark_monit_eng').switch(DesignToField).join(Field)
            for t in field: monit_fields.append(t.field.field_id)
            monit_fields=list(set(monit_fields))
            tmonit_fields=np.char.add(np.asarray(monit_fields).astype(str),'\n').tolist()
            if ptt.exists(ptt.join(daily_dir,'etc')):
                with open(ptt.join(daily_dir,'etc', 'DarkMonitor_fields'), 'w') as f:
                    f.writelines(tmonit_fields)
        except:
            monit_fields = []
            if not cron: warnings.warn('DarkMonitor_fields file is missing and DB query failed', UserWarning)
            else: splog.warnings('Warning: DarkMonitor_fields file is missing and DB query failed')
    else:
        with open(ptt.join(daily_dir,'etc', 'DarkMonitor_fields'), 'r') as f:
            monit_fields = np.array(f.readlines())
            monit_fields= np.char.replace(monit_fields, '\n','').astype(int).tolist()
    return (rm_fields, monit_fields)

def plot_milestone(obs, axs, max_mjd, im=3, jm=2, html=False):
    milestones = read_table_yanny(ptt.join(daily_dir,'etc','fiber_milestones.par'), 'MILESTONE')
    milestones = milestones[np.where(milestones['obs'] == obs)]
    labels = []
    for row in milestones:
        for i in range(0,im):
            for j in range(0,jm):

                if row['mjd'] > max_mjd + 1:
                    continue
                label = row['label'] if row['label'] != '' else None

                if not html:
                    if row['linestyle'] in ['-', '--', '-.', ':', 'None', ' ',
                                            '', 'solid', 'dashed', 'dashdot', 'dotted']:
                        ls = row['linestyle']
                    else:
                        ls = ( 0, tuple(np.array(tuple(row['linestyle'].replace(',','').replace(' ',''))).astype(int)))
                    if jm == 1:
                        axs[i].axvline(row['mjd'], color=row['color'], lw = .5, ls = ls, label=label)
                    else:
                        axs[i,j].axvline(row['mjd'], color=row['color'], lw = .5, ls = ls, label=label)
                else:
                    if row['bklinestyle'] is None: continue
                    #'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
                    bk2p = {'dotted':'dot','dashed':'dash','3 3 1 3 1 3':'longdashdot'}
                    if row['bklinestyle']  in bk2p.keys():
                        row['bklinestyle']  = bk2p[row['bklinestyle']]

                    taxs = dict(row = i+1, col=1) if jm == 1 else dict(row=i+1,col=j+1)
                    line = dict(color=colors[row['color']], width=1, dash=row['bklinestyle'])
                    axs.add_shape(type="line", x0=row['mjd'],x1=row['mjd'],y0=0,y1=1, opacity=.8,
                              xref=f'x1', yref=f'y domain', name=label, showlegend=False,
                              line=line,**taxs)
                    # Add an invisible scatter trace to serve as a legend entry for the line
                    if (label is not None and i == 0 and j == 0):
                        axs.add_trace(go.Scatter(x=[None],y=[None],mode='lines', opacity=.8,
                                      line=line,showlegend=True, name=label))
    return(axs, milestones)

def split_by_nan(x, y):
    isnan = np.isnan(y)
    segments = []
    start_idx = 0
    for idx in range(1, len(isnan)):
        if isnan[idx] != isnan[idx - 1]:
            if not isnan[start_idx:idx].all():
                segments.append((x[start_idx:idx], y[start_idx:idx]))
            start_idx = idx
    if not isnan[start_idx:].all():
        segments.append((x[start_idx:], y[start_idx:]))
    return segments
    
def plot_QA(run2ds, test, mjds={}, obs='APO', testp='/test/sean/', clobber_lists=False,
            publish = False, epoch=False, cron = False, fast_opsdb=False, html=False,
            html_name = None):

    if cron:
        makedirs(ptt.join(daily_dir,'logs','QA'),exist_ok=True)
        logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
        logging.getLogger('peewee').setLevel(logging.WARNING)
        splog = logging.getLogger('root')
        filelog = logging.FileHandler(ptt.join(daily_dir,'logs','QA',
                            datetime.datetime.today().strftime("%m%d%Y")+f'-{obs}.log'))
        filelog.setLevel(logging.DEBUG)
        filelog.setFormatter(Formatter())
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        console.setFormatter(Formatter())
        splog.addHandler(filelog)
        splog.addHandler(console)
        splog.setLevel(logging.DEBUG)
    if plotly is None and html is True:
        if cron:
            splog.warnings('plotly not installed... defaulting to no html')
        else:
            warnings.warn('plotly not installed... defaulting to no html', UserWarning)
        html = False
            
    rm_fields, monit_fields = load_fields(clobber_lists=clobber_lists)
    all_data = None
    all_mdate = ''
    all_mdate_f = 0

    for r,run2d in enumerate(run2ds):
        if test[r] is True: test_path = testp
        else: test_path = ''
        old_paths = False

        es = 'epoch' if epoch else 'daily'
        try:
            datafile = ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path, run2d,'summary',es,'spCalib_QA-'+run2d+'.fits')
            data = Table.read(datafile, format='fits')
        except:
            print('Checking old path')
            datafile = ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path, run2d,'spCalib_QA-'+run2d+'.csv')
            data=pd.read_csv(datafile)
            old_paths = True
        mdate_f = ptt.getctime(datafile)
        mdate = time.ctime(ptt.getctime(datafile))
        if all_mdate_f < mdate_f:
            all_mdate_f = mdate_f
            all_mdate = mdate
        if mjds is not None:
            if run2d in mjds.keys():
                if mjds[run2d][0] is not None: data=data[data['MJD']>=mjds[run2d][0]]
                if mjds[run2d][1] is not None: data=data[data['MJD']<=mjds[run2d][1]]
        if obs is not None:
            data=data[data['OBS'] == obs]
        max_mjd = np.max(data['MJD'].values)
        RM=data[data['FIELD'].isin(rm_fields)]
        Monit=data[data['FIELD'].isin(monit_fields)]
        if all_data is None:
            all_data = data
        else:
            all_data = pd.concat([all_data,data], ignore_index=True)
        if not html:
            rm_style    = {'ls':'', 'marker':'.', 'color':'C1', 'alpha':.5, 'label':'RM Fields'}
            main_style  = {'ls':'', 'marker':'.', 'color':'C0', 'alpha':.2, 'label':'All Fields'}
            monit_style = {'ls':'', 'marker':'.', 'color':'C3', 'alpha':.5, 'label':'DarkMonitoring Fields'}
            ylim = [-0.10,0.10]
            ylim_r = [-0.001,0.15]
            if r == 0:
                fig, axs = plt.subplots(3,2, figsize=[12,6])
                axs, milestones = plot_milestone(obs, axs, max_mjd)
            else:
                rm_style['label'] = None
                main_style['label'] = None
                monit_style['label'] = None
            for i in range(0,3):
                axs[i,0].axhline(0, alpha= .2, color='k')
                axs[i,1].axhline(0.025, alpha= .2, color='k')
                axs[i,0].plot(data['MJD'], data[filters[i]+'_MEAN'], **main_style)
                axs[i,0].plot(Monit['MJD'], Monit[filters[i]+'_MEAN'], **monit_style)
                axs[i,0].plot(RM['MJD'], RM[filters[i]+'_MEAN'], **rm_style)
                axs[i,1].plot(data['MJD'], data[filters[i]+'_SIG'], **main_style)
                axs[i,1].plot(Monit['MJD'], Monit[filters[i]+'_SIG'], **monit_style)
                axs[i,1].plot(RM['MJD'], RM[filters[i]+'_SIG'], **rm_style)
                axs[i,0].set_ylim(ylim)
                axs[i,1].set_ylim(ylim_r)
                if i != 2:
                    axs[i,0].set_xticklabels([])
                    axs[i,1].set_xticklabels([])
                else:
                    axs[i,0].set_xlabel('MJD')
                    axs[i,1].set_xlabel('MJD')
                axs[i,0].set_ylabel('['+filters[i].lower()+']')
                axs[i,1].set_ylabel('['+filters[i].lower()+']')
                if i == 0:
                    axs[0,0].set_title('mean log(synflux/calibflux)')
                    axs[0,1].set_title('sigma log(synflux/calibflux)')
        else:
            if r == 0:
                test_path = testp if test[0] is True else ''
                save_dir = getenv('BOSS_QA_DIR', default=getenv('BOSS_SPECTRO_REDUX'))
                if html_name is None:
                    html_name = f'BOSS_QA-{obs}.html'
                savename = ptt.join(save_dir,html_name)
                #output_file(filename=savename, title=f'{obs} BOSS QA')
                axs = make_subplots(rows=4, cols=2,shared_xaxes=True, shared_yaxes=False,
                                    vertical_spacing=0.02,horizontal_spacing=.05,
                                    subplot_titles = ('mean log(synflux/calibflux)',
                                                      'sigma log(synflux/calibflux)',
                                                      '','','','','',''))
                axs.update_xaxes(title='MJD',row=4)
                for r in [1,2,3]:
                    axs.update_yaxes(title=f'[{filters[r-1].lower()}]',row=r)
                axs.update_yaxes(title='SEEING50',row=4)
                axs.update_layout( paper_bgcolor='white', plot_bgcolor='white')
                for fr in [1,2,3]:
                    axs.update_yaxes(range=[-0.10,0.10],  row=fr, col=1)
                    axs.update_yaxes(range=[-0.001,0.15], row=fr, col=2)
                    axs.update_xaxes(showgrid=False, row=fr, **axopts)
                    axs.update_yaxes(showgrid=False, row=fr, **ayopts)
                    
                    axs.add_shape(type='line', x0=0, x1=1, y0=0.025, y1=0.025,
                                    xref='x domain', yref=f'y1',opacity=.2,
                                    line=dict(color='black', dash='solid',),
                                    row=fr, col=2)
                    axs.add_shape(type='line', x0=0, x1=1, y0=0, y1=0,
                                    xref='x domain', yref=f'y1',opacity=.2,
                                    line=dict(color='black', dash='solid',),
                                    row=fr, col=1)

                axs.update_xaxes(showgrid=True, row=4, **axopts)
                axs.update_yaxes(showgrid=True, row=4, **ayopts)
                axs, milestones = plot_milestone(obs, axs, max_mjd, im=4, jm=2, html=True)
            rm_style    = dict(marker=dict(size=5, color='orange'),name='RM Fields',
                               mode='markers', opacity=.5)

            main_style  = dict(marker=dict(size=5, color='blue'),name='All Fields',
                               mode='markers', opacity=.1)
                                        
            monit_style = dict(marker=dict(size=5, color='red'), mode='markers',
                               name='DarkMonitoring Fields', opacity=.5)

            for fr  in [1,2,3]:
                for c in [1,2]:
                    sl = True if (fr == 1 and c==1 and r == 0) else False

                    type = '_MEAN' if c ==1 else '_SIG'
                    ht='MJD,Mean: %{x:.2f},%{y:.2f}' if c ==1 else 'MJD,Sigma: %{x:.2f},%{y:.2f}'
                    if len(data) > 0:
                        axs.add_trace(go.Scatter(x=data['MJD'], y=data[filters[fr-1]+type],
                                      showlegend=sl, **main_style,hovertemplate=ht), row=fr,col=c)
                    if len(Monit) > 0:
                        axs.add_trace(go.Scatter(x=Monit['MJD'], y=Monit[filters[fr-1]+type],
                                      showlegend=sl, **monit_style,hovertemplate=ht), row=fr,col=c)
                    if len(RM) > 0:
                        axs.add_trace(go.Scatter(x=RM['MJD'], y=RM[filters[fr-1]+type],
                                      showlegend=sl, **rm_style,hovertemplate=ht), row=fr,col=c)

    if not html:
        if publish:
            plt.tight_layout(rect=(0,.07,1,1))
            nmiles = len(milestones[milestones['label'] != ''])
            axs[2,0].legend(frameon=True,fontsize=8,loc=2,ncol=((nmiles+3)//2), bbox_to_anchor=(0, -.35))
        else:
            plt.tight_layout()
            axs[0,0].legend(frameon=True,fontsize=3,loc=2,ncol=2)
            plt.gcf().text(0.02,0.02,'Updated: '+mdate+ ' (Latest MJD:'+str(max(data['MJD']))+')', fontsize=6, ha='left')
            plt.gcf().text(0.98,0.02,'Run2d='+','.join(run2ds), fontsize=6, ha='right')
            plt.gcf().text(0.50,0.02,'Observatory='+obs, fontsize=6, ha='center')
    
        ctype = 'epoch' if epoch else 'daily'
        if len(run2ds) == 1:
            outdir = ptt.join(getenv("BOSS_SPECTRO_REDUX"), test_path, run2d, 'summary')
            outdir = ptt.join(outdir,ctype,'spCalib_QA-'+run2ds[0]+'-'+obs+'.png')
        else:
            outdir = ptt.join(getenv("BOSS_SPECTRO_REDUX"), test_path)
            outdir = ptt.join(outdir,f'spCalib_QA-{ctype}-'+'+'.join(run2ds)+'-'+obs+'.png')

        plt.savefig(outdir, dpi=200)
        plt.show()
        if len(run2ds) == 1:
            plotsn2=True
        else:
            plotsn2 = False
    else:
        see = None
        dark_style   = dict(marker=dict(size=3, color='orange'),name='Dark',
                            mode='markers', opacity=.5)
        bright_style = dict(marker=dict(size=3, color='green'),name='Bright',
                            mode='markers', opacity=.5)
        plate_style  = dict(marker=dict(size=3, color='blue'),name='Plate',
                            mode='markers', opacity=.5)
        bidx = []
        didx = []
        data = all_data
        if (len(run2ds) == 1) or (len(mjds) == len(run2ds)):
            plotsn2 = True
        else:
            plotsn2 = False

    if plotsn2:
        mjd_limits = mjds
        spall = None
        mjds = []
        mj = []
        fsn2g = []
        fsn2r = []
        fsn2i = []
        exptime = []
        fcad = []
        see = []
        idx = 0

        ep = '' if not epoch else 'epoch'
        if cron: splog.info(f"Reading Field-MJD spAll for run2d={' '.join(run2ds)}")

        for i, row in tqdm(data.iterrows(), total=data.shape[0]):
            idx = idx + 1
            if mjd_limits is not None:
                r2d = None
                for ir2d, run2d in enumerate(run2ds):
                    if mjd_limits[run2d][0] is not None:
                        if int(row['MJD']) < mjd_limits[run2d][0]:
                            continue
                    if mjd_limits[run2d][1] is not None:
                        if int(row['MJD']) > mjd_limits[run2d][1]:
                            continue
                    r2d = run2d
                    break
            else:
                ir2d = 0
                r2d = run2ds[0]
            if r2d is None:
                if cron:
                    splog.info('skipping: spAll-'+str(row['FIELD']).zfill(6)+'-'+str(row['MJD'])+'.fits')
                continue

            test_path = testp if test[ir2d] is True else ''
            if not old_paths:
                spallfile = glob(ptt.join(field_spec_dir(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path),
                                                run2d, row['FIELD'], row['MJD'],
                                                epoch = epoch, full= True),
                                                f"spAll-{str(row['FIELD']).zfill(6)}-{row['MJD']}.fits*"))
            else:
                if epoch:
                    spallfile = glob(ptt.join(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path),
                                                    run2d, 'epoch','spectra','full',
                                                    str(row['FIELD']).zfill(6), str(row['MJD']),
                                                   f"spAll-{str(row['FIELD']).zfill(6)}-{row['MJD']}.fits*"))

                else:
                    spallfile = glob(ptt.join(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path),
                                                    run2d, 'spectra','full',
                                                    str(row['FIELD']).zfill(6), str(row['MJD']),
                                                   f"spAll-{str(row['FIELD']).zfill(6)}-{row['MJD']}.fits*"))

            if cron:
                splog.info(f'({idx}/{data.shape[0]}) '+'spAll-'+str(row['FIELD']).zfill(6)+'-'+str(row['MJD'])+'.fits'+ f'(run2d:{r2d})')
            if len(spallfile) == 0:
                continue
            else:
                spallfile = spallfile[0]
            try:
                spall = Table(fits.getdata(spallfile,1))
            except Exception as e:
                print(e)
                continue
            for row in (spall):
                tobs = None
                if np.isnan(row['MJD_FINAL']):
                    continue
                if row['MJD_FINAL'] <= 59847:
                    tobs = 'APO'
                elif (int(row['FIELD']) > 100000) & (row['MJD_FINAL'] < 60191):
                    tobs = 'APO'
                else:
                    if 'OBS' in row.colnames:
                        tobs = row['OBS']
                    else: 
                        try:
                            tobs = row['OBSERVATORY']
                        except:
                           continue 
                if tobs not in ['APO','LCO']:
                    idx = np.where(data['FIELD'] == int(row['FIELD']))[0]
                    if len(idx) == 0: continue
                    else: tobs = obs.upper()
                if tobs != obs.upper():
                    continue
                nexp = len(row['FIELDSNR2R_LIST'].split())

                if row['MJD_FINAL'] == 0:
                     mjds.extend([row['MJD']]*nexp)
                else:
                     mjds.extend([row['MJD_FINAL']]*nexp)

#                mj.extend([row['MJD']]*nexp)
                mj.extend((np.atleast_1d(np.asarray(row['TAI_LIST'].split())).astype(float)/(24*3600)).tolist())
                fsn2g.extend(row['FIELDSNR2G_LIST'].split())
                fsn2r.extend(row['FIELDSNR2R_LIST'].split())
                fsn2i.extend(row['FIELDSNR2I_LIST'].split())
                exptime.extend([row['EXPTIME']/row['NEXP']]*nexp)
                fcad.extend([row['CADENCE']]*nexp)
                see.extend([row['SEEING50']]*nexp)
        mjds = np.asarray(mjds)
        mj   = np.asarray(mj)
        exptime = np.asarray(exptime)
        fcad = np.asarray(fcad)
        fsn2g = np.asarray(fsn2g).astype(float)
        fsn2r = np.asarray(fsn2r).astype(float)
        fsn2i = np.asarray(fsn2i).astype(float)
        see = np.asarray(see).astype(float)
        test = Table({'MJD': mjds,'MJ':mj, 'exptime':exptime, 'fcad':fcad,'fsn2g':fsn2g,'fsn2r':fsn2r,'fsn2i':fsn2i,'see':see})
        if len(test) == 0:
            return
        test = unique(test,keys=['MJ','exptime','fcad','fsn2g','fsn2r','fsn2i'])
        mjds = test['MJ'].data
        exptime = test['exptime'].data
        fsn2g = test['fsn2g'].data
        fsn2r = test['fsn2r'].data
        fsn2i = test['fsn2i'].data
        fcad = test['fcad'].data
        see = test['see'].data

        didx = np.where(wwhere(fcad, 'dark*'))[0]
        bidx = np.where(wwhere(fcad, 'bright*'))[0]
        pidx = np.where(wwhere(fcad, 'plate*'))[0]
        
        if not html:
            fig, axs = plt.subplots(4,1, figsize=[12,8])

            axs, milestones = plot_milestone(obs, axs, max_mjd, im=4, jm=1)
        
            axs[0] = plot_sn2_filt(axs[0], mjds, fsn2g, exptime, pidx, bidx, didx, 'FieldSN2_G/(900s)')
            axs[1] = plot_sn2_filt(axs[1], mjds, fsn2r, exptime, pidx, bidx, didx, 'FieldSN2_R/(900s)')
            axs[1].set_ylim(-.5,22)
            axs[2] = plot_sn2_filt(axs[2], mjds, fsn2i, exptime, pidx, bidx, didx, 'FieldSN2_I/(900s)')
            axs[2].set_ylim(-.5,22)
            
            
            axs[3] = plot_sn2_filt(axs[3], mjds, see, np.full_like(see, 900.0), pidx, bidx, didx, 'SEEING50', labelbottom=True)
            ylim = list(axs[3].get_ylim())
            if ylim[1] > 5.2:
                ylim[1] = 5.2
            if ylim[0] > 0:
                ylim[0] = -.2
            axs[3].set_ylim(ylim)

            axs[0].set_title('Field SNR^2 Per Exposure (900s)')
            fig.tight_layout()
            plt.subplots_adjust(hspace=0,bottom=.08)
            plt.gcf().text(0.02,0.02,'Updated: '+mdate+ ' (Latest MJD:'+str(max(mj))+')', fontsize=6, ha='left')
            plt.gcf().text(0.98,0.02,'Run2d='+','.join(run2ds), fontsize=6, ha='right')
            plt.gcf().text(0.50,0.02,'Observatory='+obs, fontsize=6, ha='center')
            outdir = ptt.join(getenv("BOSS_SPECTRO_REDUX"), test_path, run2d)
            outdir = ptt.join(outdir,'summary')
            if not epoch:
                plt.savefig(ptt.join(outdir,'daily','SN2-'+run2ds[0]+'-'+obs+'.png'), dpi=200)
            else:
                plt.savefig(ptt.join(outdir,'epoch','SN2-'+run2ds[0]+'-'+obs+'.png'), dpi=200)
            plt.show()
        else:
            axs2 = make_subplots(rows=4, cols=1, shared_xaxes=True, shared_yaxes=False,vertical_spacing=0.02,
                                 subplot_titles=('Field SNR^2 Per Exposure (900s)','','',''))
            axs2.update_xaxes(title='MJD',row=4)
            for r in [1,2,3]:
                axs2.update_yaxes(title=f'FieldSN2_{filters[r-1].upper()}/(900s)',row=r)
                axs2.update_xaxes(showgrid=True, row=r, **axopts)
                axs2.update_yaxes(showgrid=True, row=r, **ayopts)
            axs2.update_xaxes(showgrid=True, row=4, **axopts)
            axs2.update_yaxes(showgrid=True, row=4, **ayopts)
            axs2.update_yaxes(title='SEEING50',row=4)
            axs2.update_layout( paper_bgcolor='white', plot_bgcolor='white')
            for r in [1,2,3]:
                axs2.add_shape(go.layout.Shape(type="line", x0=0,x1=1,y0=0,y1=0,
                                              xref=f'x domain', yref=f'y', opacity=.2,
                                              line=dict(color="black", dash="solid")),
                                              row=r,col=1)

            axs2, milestones = plot_milestone(obs, axs2, max_mjd, im=4, jm=1, html=True)
            fsn2 = {'G':fsn2g,'R':fsn2r,'I':fsn2i}
            didx = np.where(wwhere(fcad, 'dark*'))[0]
            bidx = np.where(wwhere(fcad, 'bright*'))[0]
            pidx = np.where(wwhere(fcad, 'plate*'))[0]
            for i, filt in enumerate(filters):
                fi = i+1
                sl = True if fi == 1 else False
                ht='MJD,FieldSN2: %{x:.2f},%{y:.2f}'
                if len(pidx) > 0:
                    axs2.add_trace(go.Scatter(x=mjds[pidx], y=fsn2[filt][pidx],
                                   showlegend=sl, **plate_style,hovertemplate=ht), row=fi, col=1)
                if len(bidx) > 0:
                    axs2.add_trace(go.Scatter(x=mjds[bidx], y=fsn2[filt][bidx],
                                   showlegend=sl, **bright_style,hovertemplate=ht), row=fi, col=1)
                if len(didx) > 0:
                    axs2.add_trace(go.Scatter(x=mjds[didx], y=fsn2[filt][didx],
                                   showlegend=sl, **dark_style,hovertemplate=ht), row=fi, col=1)
                if len(didx) > 0:
                    moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],fsn2[filt][didx]*900.0/exptime[didx])
                    color=dark_style['marker']['color']
                    for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                        axs2.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                            line=dict(color=color, width=.5), showlegend=False),row=fi, col=1)
                        axs2.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                            mode='lines', line=dict(color=color, width=.5), showlegend=False), row=fi, col=1)
                    axs2.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                            line=dict(color=color, width=1),showlegend=False), row=fi, col=1)
                if len(pidx) > 0:
                    moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],fsn2[filt][pidx]*900.0/exptime[pidx])
                    color=plate_style['marker']['color']
                    for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                        axs2.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                            line=dict(color=color, width=.5), showlegend=False),row=fi, col=1)
                        axs2.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                            mode='lines', line=dict(color=color, width=.5), showlegend=False), row=fi, col=1)
                    axs2.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                            line=dict(color=color, width=1),showlegend=False), row=fi, col=1)

            i = 3
            fr = i+1
            ht='MJD,SEEING50: %{x:.2f},%{y:.2f}'
            if len(pidx) > 0:
                axs2.add_trace(go.Scatter(x=mjds[pidx], y=see[pidx],
                               showlegend=False, **plate_style,hovertemplate=ht), row=fr, col=1)
                for c in [1,2]:
                    axs.add_trace(go.Scatter(x=mjds[pidx], y=see[pidx],showlegend=False,
                                  **plate_style,hovertemplate=ht), row=4, col=c)
            if len(bidx) > 0:
                axs2.add_trace(go.Scatter(x=mjds[bidx], y=see[bidx],
                               showlegend=False, **bright_style,hovertemplate=ht), row=fr, col=1)
                for c in [1,2]:
                    axs.add_trace(go.Scatter(x=mjds[bidx], y=see[bidx],showlegend=False,
                                  **bright_style,hovertemplate=ht), row=4, col=c)
            if len(didx) > 0:
                axs2.add_trace(go.Scatter(x=mjds[didx], y=see[didx],
                               showlegend=False, **dark_style,hovertemplate=ht), row=fr, col=1)
                for c in [1,2]:
                    axs.add_trace(go.Scatter(x=mjds[didx], y=see[didx],showlegend=False,
                                  **dark_style,hovertemplate=ht), row=4, col=c)

            if len(didx) > 0:
                moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],see[didx])
                color=dark_style['marker']['color']
                for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                    axs2.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                        line=dict(color=color, width=.5), showlegend=False),row=4, col=1)
                    axs2.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                        mode='lines', line=dict(color=color, width=.5), showlegend=False), row=4, col=1)
                axs2.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                        line=dict(color=color, width=1),showlegend=False), row=4, col=1)


                for c in [1,2]:
                    for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                        axs.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                            line=dict(color=color, width=.5), showlegend=False),row=4, col=c)
                        axs.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                            mode='lines', line=dict(color=color, width=.5), showlegend=False), row=4, col=c)
                    axs.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                            line=dict(color=color, width=1),showlegend=False), row=4, col=c)

            if len(pidx) > 0:
                moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],see[pidx])
                color=plate_style['marker']['color']
                for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                    axs2.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                        line=dict(color=color, width=.5), showlegend=False),row=4, col=1)
                    axs2.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                        mode='lines', line=dict(color=color, width=.5), showlegend=False), row=4, col=1)
                axs2.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                        line=dict(color=color, width=1),showlegend=False), row=4, col=1)

                for c in [1,2]:
                    for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                        axs.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                            line=dict(color=color, width=.5), showlegend=False),row=4, col=c)
                        axs.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                            mode='lines', line=dict(color=color, width=.5), showlegend=False), row=4, col=c)
                    axs.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                            line=dict(color=color, width=1),showlegend=False), row=4, col=c)

            axs.update_layout(legend=dict(orientation='h', x=0.5,y=-0.1, xanchor='center',
                                         yanchor='top', font=dict(size=8),traceorder='normal',
                                         itemwidth=30, tracegroupgap=1,borderwidth=.5,
                                         bgcolor='rgba(255, 255, 255, 0.5)', bordercolor='black'))
            axs2.update_layout(legend=dict(orientation='h', x=0.5,y=-0.12, xanchor='center',
                                         yanchor='top', font=dict(size=8),traceorder='normal',
                                         itemwidth=30, tracegroupgap=1,borderwidth=.5,
                                         bgcolor='rgba(255, 255, 255, 0.5)', bordercolor='black'))

#        axs.update_layout(legend=dict(x=0,y=1,traceorder='normal',orientation='v',
#                                     font=dict(size=8), bgcolor='rgba(255, 255, 255, 0.5)',
#                                     bordercolor='black', borderwidth=.5, tracegroupgap=1))
#        axs2.update_layout(legend=dict(x=0,y=1,traceorder='normal',orientation='v',
#                                     font=dict(size=8), bgcolor='rgba(255, 255, 255, 0.5)',
#                                     bordercolor='black', borderwidth=.5, tracegroupgap=1))
            if see is not None:
                ymin = np.nanmin(see)
                ymax = np.nanmax(see)
                if ymin < 0:
                    ymin = -0.2
                if ymax > 5.2:
                    ymax = 5.2
                axs.update_yaxes(range=[ymin, ymax], row=4)
                axs2.update_yaxes(range=[ymin, ymax], row=4)


            if cron: splog.info(f'Reading OPSDB Exposures')
            try:
                opsdb.database.set_profile('operations')
                environ['OBSERVATORY'] = obs.upper()
                opsdb.database.connect()  # This will recreate the base model class and reload all the model classes
                CCDs = ['b2','r2'] if obs.upper() == 'LCO' else ['b1','r1']
                from sdssdb.peewee.sdss5db.opsdb import Exposure, CameraFrame, Camera, Configuration
            except Exception as e:
                print(e)
                if not cron: warnings.warn('No Connection to the SDSS-V OpsDB (internal only)', UserWarning)
                else: splog.warnings('WARNING: No Connection to the SDSS-V OpsDB (internal only)')
                sos_data = None
                Configuration = None
            
            if Configuration is not None:
                if not fast_opsdb:
                    blue = CameraFrame.select(CameraFrame.exposure.exposure_no.alias('expid'),CameraFrame.exposure.exposure_time.alias('exptime'), CameraFrame.exposure.start_time, CameraFrame.sn2.alias(CCDs[0]), Design.design_mode.alias('design_mode')).join(Exposure).join(Configuration).join(Design,  on=(Configuration.design_id == Design.design_id)).join(DesignMode).switch(CameraFrame).join(Camera).where(Camera.label == CCDs[0]).dicts()
                    red  = CameraFrame.select(CameraFrame.exposure.exposure_no.alias('expid'), CameraFrame.sn2.alias(CCDs[1])).join(Exposure).switch(CameraFrame).join(Camera).where(Camera.label == CCDs[1]).dicts()
                    sos_data=pd.DataFrame(blue)
                    sos_data_r = pd.DataFrame(red)
                    sos_data = pd.merge(sos_data, sos_data_r, on='expid', how='outer', indicator=True)
                    sos_data = sos_data.loc[sos_data['_merge'] == 'both']
                    sos_data = sos_data.assign(mjd=lambda x: Time(x.start_time).mjd)
                    sos_data = sos_data.sort_values('mjd')
                    _dir = ptt.dirname(savename)
                    _name = ptt.basename(savename).replace('.html','_sos.csv')
                    sos_data.to_csv(ptt.join(_dir,f'.{_name}'))
                else:
                    _dir = ptt.dirname(savename)
                    _name = ptt.basename(savename).replace('.html','_sos.csv')
                    if ptt.exists(ptt.join(_dir,f'.{_name}')):
                        sos_data = pd.read_csv(ptt.join(_dir,f'.{_name}'))
                    else:
                        sos_data = None
                        if cron:splog.info(f'No OPSDB SOS data')
                        else: warnings.warn('No OPSDB SOS data')
                if cron: splog.info(f'Done Reading OPSDB Exposures')
            
            if sos_data is not None:
                for ccd in CCDs:
                    sos_data[f'{ccd}_900'] = sos_data[ccd]/sos_data['exptime']*900
                axs3 = make_subplots(rows=3, cols=1, shared_xaxes=True, shared_yaxes=False,vertical_spacing=0.02,
                                     subplot_titles=('SOS SNR^2 Per Exposure (900s)','','',''))
                axs3.update_xaxes(title='MJD',row=3)
                axs3.update_yaxes(title=f'SOS SN2_{CCDs[0]}/(900s)',row=1)
                axs3.update_yaxes(title=f'SOS SN2_{CCDs[1]}/(900s)',row=2)
                axs3.update_yaxes(title='SEEING50',row=3)
                axs3.update_layout( paper_bgcolor='white', plot_bgcolor='white')
                for r in [1,2,3]:
                    axs3.update_xaxes(showgrid=True, row=r, **axopts)
                    axs3.update_yaxes(showgrid=True, row=r, **ayopts)
                axs3, milestones = plot_milestone(obs, axs3, max_mjd, im=3, jm=1, html=True)
                dsidx = np.where(wwhere(sos_data['design_mode'].values, 'dark*'))[0]
                bsidx = np.where(wwhere(sos_data['design_mode'].values, 'bright*'))[0]

                for i, ccd in enumerate(CCDs):
                    r = i + 1
                    sl = True if r == 1 else False
                    ht = 'MJD,SOS_SN2: %{x:.2f},%{y:.2f}'
                    axs3.add_shape(go.layout.Shape(type="line", x0=0,x1=1,y0=0,y1=0,
                                                  xref=f'x domain', yref=f'y', opacity=.2,
                                                  line=dict(color="black", dash="solid")),
                                                  row=r,col=1)
                    
                    if len(bsidx) > 0:
                        axs3.add_trace(go.Scatter(x=sos_data['mjd'].values[bsidx],
                                                  y=sos_data[f'{ccd}_900'].values[bsidx],
                                       showlegend=sl, **bright_style,hovertemplate=ht), row=r, col=1)

                    if len(dsidx) > 0:
                        axs3.add_trace(go.Scatter(x=sos_data['mjd'].values[dsidx],
                                                  y=sos_data[f'{ccd}_900'].values[dsidx],
                                       showlegend=sl, **dark_style,hovertemplate=ht), row=r, col=1)
                        color=dark_style['marker']['color']
                        moving_mjd, moving_avg, moving_16, moving_84 = runAvg(sos_data['mjd'].values[dsidx],
                                                                              sos_data[f'{ccd}_900'].values[dsidx])
                        
                        for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                            axs3.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                                line=dict(color=color, width=.5), showlegend=False),row=r, col=1)
                            axs3.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                                mode='lines', line=dict(color=color, width=.5), showlegend=False), row=r, col=1)
                        axs3.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                                line=dict(color=color, width=1),showlegend=False), row=r, col=1)

                i = 2
                fr = i+1
                ht = 'MJD,Seeing: %{x:.2f},%{y:.2f}'
                if len(bidx) > 0:
                    axs3.add_trace(go.Scatter(x=mjds[bidx], y=see[bidx],
                                   showlegend=False, **bright_style,hovertemplate=ht), row=fr, col=1)
                if len(didx) > 0:
                    axs3.add_trace(go.Scatter(x=mjds[didx], y=see[didx],
                                   showlegend=False, **dark_style,hovertemplate=ht), row=fr, col=1)
                    moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],see[didx])
                    color=dark_style['marker']['color']
                    for (x_seg, upper_seg), (_, lower_seg) in zip(split_by_nan(moving_mjd, moving_84), split_by_nan(moving_mjd, moving_16)):
                        axs3.add_trace(go.Scatter(x=x_seg, y=upper_seg,fill=None, mode='lines',opacity=.1,
                            line=dict(color=color, width=.5), showlegend=False),row=fr, col=1)
                        axs3.add_trace(go.Scatter(x=x_seg, y=lower_seg, fill='tonexty',opacity=.1,
                            mode='lines', line=dict(color=color, width=.5), showlegend=False), row=fr, col=1)
                    axs3.add_trace(go.Scatter(x=moving_mjd,y=moving_avg, mode='lines',opacity=1,
                            line=dict(color=color, width=1),showlegend=False), row=fr, col=1)

    #            axs3.update_layout(legend=dict(x=0,y=1,traceorder='normal',orientation='v',
    #                                         font=dict(size=8), bgcolor='rgba(255, 255, 255, 0.5)',
    #                                         bordercolor='black', borderwidth=.5, tracegroupgap=1))
                axs3.update_layout(legend=dict(orientation='h', x=0.5,y=-0.12, xanchor='center',
                                             yanchor='top', font=dict(size=8),traceorder='normal',
                                             itemwidth=30, tracegroupgap=1,borderwidth=.5,
                                             bgcolor='rgba(255, 255, 255, 0.5)', bordercolor='black'))

                if see is not None:
                    ymin = np.nanmin(see)
                    ymax = np.nanmax(see)
                    if ymin < 0:
                        ymin = -0.2
                    if ymax > 5.2:
                        ymax = 5.2
                    axs3.update_yaxes(range=[ymin, ymax], row=4)
                sosmjd = max(sos_data['mjd'])
            else:
                axs3= []
                sosmjd = ''
            
            title = '\n'.join([f"<h2>SDSS-V {obs} BOSS QA: {','.join(run2ds)}</h2>",
                               f"<p>Spectro-photometric Plot Updated: {all_mdate}<br>",
                               f"   Latest Full Pipeline MJD:{max(data['MJD'])}<br>",
                               f"   Latest SOS Pipeline MJD:{sosmjd}<br>",
                               f"   Last Updated: {time.ctime()} (MJD: {int(Time.now().mjd)})</p>"])

            config1 = {'toImageButtonOptions': {'format': 'png','scale': 6 , # Multiply title/legend/axis/canvas sizes by this factor
                                               'filename': f"spCalib_QA-{','.join(run2ds)}-{obs}.png"},
                      'responsive': True}  # Ensure the figure is responsive
            fig1_params = dict(full_html=False, default_height='900px', default_width='100%',
                              include_plotlyjs='cdn', config=config1)

            config2 = {'toImageButtonOptions': {'format': 'png','scale': 6 , # Multiply title/legend/axis/canvas sizes by this factor
                                               'filename': f"FieldSN2-{','.join(run2ds)}-{obs}.png"},
                      'responsive': True}  # Ensure the figure is responsive
            fig2_params = dict(full_html=False, default_height='900px', default_width='100%',
                              include_plotlyjs='cdn', config=config2)
            plotly_jinja_data = {"fig": axs.to_html(**fig1_params),
                                 "fig_FieldSN2": axs2.to_html(**fig2_params),
                                 "title":title, "name":f'{obs} BOSS QA', "favicon":favicon}

            if see is not None:
                config3 = {'toImageButtonOptions': {'format': 'png','scale': 6 , # Multiply title/legend/axis/canvas sizes by this factor
                                                   'filename': f"SOS_SN2-{','.join(run2ds)}-{obs}.png"},
                           'responsive': True}  # Ensure the figure is responsive
                fig3_params = dict(full_html=False, default_height='900px', default_width='100%',
                                  include_plotlyjs='cdn', config=config3)

                plotly_jinja_data["fig_SOS"] = axs3.to_html(**fig3_params)
                template = ptt.join(idlspec2d_dir,'templates','html','QA_template.html')
            else:
                template = ptt.join(idlspec2d_dir,'templates','html','QA_noSOS_template.html')

            with open(savename, "w", encoding="utf-8") as output_file:
                with open(template) as template_file:
                    j2_template = Template(template_file.read())
                    output_file.write(j2_template.render(plotly_jinja_data))
                    

def plot_sn2_filt(axs, mjds, fsn2, exptime, pidx,bidx,didx,label, labelbottom=False):

    axs.grid(linestyle = ':', linewidth = 0.5)
    axs.minorticks_on()
    axs.tick_params(axis='both', which='both', labelleft=True, labelright=True, labelbottom=labelbottom, labeltop=False, bottom=True, top=False, left=True, right=True)
    axs.tick_params(axis='x', which='both', bottom=True, top=True, direction='inout')
    axs.plot(mjds[pidx], fsn2[pidx]*900.0/exptime[pidx], ls = '', marker= '.', ms = .5, color='C0', alpha=.5)
    axs.plot(mjds[bidx], fsn2[bidx]*900.0/exptime[bidx], ls = '', marker= '.', ms = .5, color='C2', alpha=.5)
    axs.plot(mjds[didx], fsn2[didx]*900.0/exptime[didx], ls = '', marker= '.', ms = .5, color='C1', alpha=.5)
    if len(didx) > 0:
        moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],fsn2[didx]*900.0/exptime[didx])
        axs.plot(moving_mjd, moving_avg, color='C1', alpha=.5)
        axs.fill_between(moving_mjd,moving_16, moving_84, color='C1', alpha=.1)
    if len(pidx) > 0:
        moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],fsn2[pidx]*900.0/exptime[pidx])
        axs.plot(moving_mjd, moving_avg, color='C0', alpha=.5)
        axs.fill_between(moving_mjd,moving_16, moving_84, color='C0', alpha=.1)
    axs.plot(np.nan, np.nan,label='plate',  ls = '', marker= '.', color='C0', ms = 2)
    axs.plot(np.nan, np.nan,label='bright', ls = '', marker= '.', color='C2', ms = 2)
    axs.plot(np.nan, np.nan,label='dark',   ls = '', marker= '.', color='C1', ms = 2)
    axs.legend(frameon=True,fontsize=3,loc=2)
    axs.set_xlabel('MJD')
    axs.set_ylabel(label)
    return(axs)

def runAvg(mjd, val, ws=7):
    i = int(np.nanmin(mjd))
    moving_mjd = []
    moving_avg = []
    moving_16  = []
    moving_84  = []
    maxi = int(np.nanmax(mjd))#-ws+1
    warnings.filterwarnings("error")

    while i < maxi:
        idx = np.where((mjd>=i) & (mjd < i+ws))[0]
        i=i+1
        if len(idx) == 0:
            moving_mjd.append(i+.5*ws)
            moving_avg.append(np.NaN)
            moving_16.append(np.NaN)
            moving_84.append(np.NaN)
        else:
            try:
                np.nanmean(mjd[idx])
                np.nanmean(val[idx])
            except:
                moving_mjd.append(i+.5*ws)
                moving_avg.append(np.NaN)
                moving_16.append(np.NaN)
                moving_84.append(np.NaN)
                continue
            moving_mjd.append(np.nanmean(mjd[idx]))
            moving_avg.append(np.nanmean(val[idx]))
            moving_16.append(np.percentile(val[idx],16))
            moving_84.append(np.percentile(val[idx],84))

    warnings.resetwarnings()

    return(moving_mjd, moving_avg, moving_16, moving_84)


