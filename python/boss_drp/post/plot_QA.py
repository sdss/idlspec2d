#!/usr/bin/env python3
from boss_drp.field import field_spec_dir
from boss_drp.utils import load_env
from boss_drp import daily_dir
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
from os import getenv, environ
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
    from bokeh.plotting import figure,output_file, save
    from bokeh.models import Div, Span, Band, ColumnDataSource, HoverTool
    try:
        from bokeh.models.widgets import Tabs, Panel
    except:
        from bokeh.models import Tabs, TabPanel as Panel
    from bokeh.palettes import viridis as colors
    from bokeh.layouts import layout, Spacer
    import bokeh
    import bs4
    colors={'g':'green','y':'yellow','r':'red','m':'magenta','b':'blue', 'k':'black'}
    tools = "box_zoom,wheel_zoom,xwheel_zoom,pan,xpan,reset,undo,redo,save"

except:
    bokeh = None



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
                    
                    taxs = axs[i] if jm == 1 else axs[i,j]
                    ls = row['bklinestyle']
                    if ls =='None': continue
                    vline = Span(location=row['mjd'], dimension='height',
                                 line_color=colors[row['color']], line_width=.5,
                                 line_dash=ls)
                    if label is not None:
                        taxs.line([np.nan], [np.nan], line_color=colors[row['color']],
                                      line_width=.5, line_dash=ls,legend_label=label)
                    if i!=0 or j!=0: taxs.legend.visible = False
                    else:
                        taxs.legend.location = "top_left"
                        taxs.legend.label_text_font_size = '4pt'
                        taxs.legend.label_height = 3
                        taxs.legend.glyph_height = 3
                        taxs.legend.spacing = 0
                        taxs.legend.background_fill_alpha = .5
                    taxs.renderers.extend([vline])

    return(axs, milestones)


def plot_QA(run2ds, test, mjds={}, obs='APO', testp='/test/sean/', clobber_lists=False,
            publish = False, epoch=False, cron = False, fast_opsdb=False, html=False,
            html_name = None):

    if cron:
        makedirs(ptt.join(daily_dir,'logs','QA'),exist_ok=True)
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
    if bokeh is None and html is True:
        if cron:
            splog.warnings('Bokeh not installed... defaulting to no html')
        else:
            warnings.warn('Bokeh not installed... defaulting to no html', UserWarning)
        html = False
            
    rm_fields, monit_fields = load_fields(clobber_lists=clobber_lists)
    all_data = None
    all_mdate = ''
    all_mdate_f = 0

    for r,run2d in enumerate(run2ds):
        if test[r] is True: test_path = testp
        else: test_path = ''
        
        es = 'epoch' if epoch else 'daily'
        csvfile = ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path, run2d,'summary',es,'spCalib_QA-'+run2d+'.csv')
        
        data=pd.read_csv(csvfile)
        mdate_f = ptt.getctime(csvfile)
        mdate = time.ctime(ptt.getctime(csvfile))
        if all_mdate_f < mdate_f:
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
                output_file(filename=savename, title=f'{obs} BOSS QA')
                axs = np.zeros([4,2],dtype=object)

                for i in range(3):
                    for j in range(2):
                        ylim = [-0.10,0.10] if j == 0 else [-0.001,0.15]
                        xlab = 'MJD' if i == 3 else None
                        if i == 0:
                            title = 'mean log(synflux/calibflux)' if j == 0 else 'sigma log(synflux/calibflux)'
                        else:
                            title = None
                        
                        if i ==0 and j == 0:
                            axs[i,j] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=f'[{filters[i].lower()}]',
                                              y_range=ylim, title=title, sizing_mode='stretch_both')
                        else:
                            axs[i,j] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=f'[{filters[i].lower()}]',
                                              y_range=ylim, title=title,x_range =axs[0,0].x_range, sizing_mode='stretch_both')
                        axs[i,j].xgrid.grid_line_color = None
                        axs[i,j].ygrid.grid_line_color = None
                        
                        tt1 = 'mean' if j == 0 else 'sigma'
                        axs[i,j].add_tools(HoverTool(tooltips=[(f'MJD,{tt1}', '$x{int},$y{0.000}')]))
                for k in [0,1]:
                    try:
                        axs[3,k] = figure(tools=tools,active_drag="box_zoom", x_axis_label='MJD',y_axis_label=f'SEEING50', title=None, x_range = axs[0,0].x_range, plot_width = axs[k,0].plot_width)
                    except:
                        axs[3,k] = figure(tools=tools,active_drag="box_zoom", x_axis_label='MJD',y_axis_label=f'SEEING50', title=None, x_range = axs[0,0].x_range, width = axs[k,0].width,sizing_mode='stretch_both')
                    axs[3,k].add_tools(HoverTool(tooltips=[('MJD,Seeing', '$x{int},$y{0.000}')]))
                    for i in range(4):
                        axs[i,k].min_border_left = 75
            axs, milestones = plot_milestone(obs, axs, max_mjd, im=4, jm=2, html=True)
            rm_style    = {'size':5, 'color':'orange', 'alpha':.5, 'legend_label':'RM Fields'}
            main_style  = {'size':5, 'color':'blue', 'alpha':.1, 'legend_label':'All Fields'}
            monit_style = {'size':5, 'color':'red', 'alpha':.5, 'legend_label':'DarkMonitoring Fields'}

            for i in range(0,3):
                hline = Span(location=0, dimension='width', line_color='black', line_alpha=.2)
                axs[i,0].renderers.extend([hline])
                hline = Span(location=0.025, dimension='width', line_color='black', line_alpha=.2)
                axs[i,1].renderers.extend([hline])
                try:
                    if len(data) > 0:  axs[i,0].scatter(data['MJD'], data[filters[i]+'_MEAN'],  **main_style)
                    if len(Monit) > 0: axs[i,0].scatter(Monit['MJD'], Monit[filters[i]+'_MEAN'], **monit_style)
                    if len(RM) > 0:    axs[i,0].scatter(RM['MJD'], RM[filters[i]+'_MEAN'],    **rm_style)
                    if len(data) > 0:  axs[i,1].scatter(data['MJD'], data[filters[i]+'_SIG'],   **main_style)
                    if len(Monit) > 0: axs[i,1].scatter(Monit['MJD'], Monit[filters[i]+'_SIG'],  **monit_style)
                    if len(RM) > 0:    axs[i,1].scatter(RM['MJD'], RM[filters[i]+'_SIG'],     **rm_style)
                except:
                    if len(data) > 0:  axs[i,0].circle(data['MJD'], data[filters[i]+'_MEAN'],  **main_style)
                    if len(Monit) > 0: axs[i,0].circle(Monit['MJD'], Monit[filters[i]+'_MEAN'], **monit_style)
                    if len(RM) > 0:    axs[i,0].circle(RM['MJD'], RM[filters[i]+'_MEAN'],    **rm_style)
                    if len(data) > 0:  axs[i,1].circle(data['MJD'], data[filters[i]+'_SIG'],   **main_style)
                    if len(Monit) > 0: axs[i,1].circle(Monit['MJD'], Monit[filters[i]+'_SIG'],  **monit_style)
                    if len(RM) > 0:    axs[i,1].circle(RM['MJD'], RM[filters[i]+'_SIG'],     **rm_style)
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
        dark_style   = {'size':1, 'color':'orange', 'alpha':.5, 'legend_label':'Dark'}
        bright_style = {'size':1, 'color':'green',  'alpha':.5, 'legend_label':'Bright'}
        plate_style  = {'size':1, 'color':'blue',   'alpha':.5, 'legend_label':'Plate'}
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
            spallfile = glob(ptt.join(field_spec_dir(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path),
                                            run2d, row['FIELD'], row['MJD'],
                                            epoch = epoch, full= True),f"spAll-{row['FIELD']}-{row['MJD']}.fits*"))
            if cron:
                splog.info(f'({idx}/{ldat}) '+'spAll-'+str(row['FIELD']).zfill(6)+'-'+str(row['MJD'])+'.fits'+ f'(run2d:{r2d})')
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

                mj.extend([row['MJD']]*nexp)
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
        mjds = test['MJD'].data
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
            axs[1] = plot_sn2_filt(axs[0], mjds, fsn2r, exptime, pidx, bidx, didx, 'FieldSN2_R/(900s)')
            axs[1].set_ylim(-.5,22)
            axs[2] = plot_sn2_filt(axs[0], mjds, fsn2i, exptime, pidx, bidx, didx, 'FieldSN2_I/(900s)')
            axs[2].set_ylim(-.5,22)
            
            
            axs[3] = plot_sn2_filt(axs[0], mjds, see, np.full_like(see, 900.0), pidx, bidx, didx, 'SEEING50')
            # np.full_like(see, 900.0) is set to create a scale factor of 1 for the seeing, as seeing does not scale with expTime
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
            axs2 = np.zeros([4],dtype=object)
            for i in range(4):
                ylim = None
                xlab = 'MJD' if i == 3 else None
                title = 'Field SNR^2 Per Exposure (900s)' if i == 0 else None
                ylab = f'FieldSN2_{filters[i]}/(900s)' if i != 3 else 'SEEING50'
                if i == 0:
                    try:
                        axs2[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab, y_range=ylim, title=title)
                    except:
                        axs2[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab, title=title,sizing_mode='stretch_both')
                else:
                    try:
                        axs2[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab, y_range=ylim, title=title, x_range = axs2[0].x_range)
                    except:
                        axs2[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab, title=title,sizing_mode='stretch_both', x_range = axs2[0].x_range)
                tt1 = 'FieldSN2' if i != 3 else 'Seeing'
                axs2[i].add_tools(HoverTool(tooltips=[(f'MJD,{tt1}', '$x{int},$y{0.000}')]))
                axs2[i].min_border_left = 75
            axs2, milestones = plot_milestone(obs, axs2, max_mjd, im=4, jm=1, html=True)
            fsn2 = {'G':fsn2g,'R':fsn2r,'I':fsn2i}
            didx = np.where(wwhere(fcad, 'dark*'))[0]
            bidx = np.where(wwhere(fcad, 'bright*'))[0]
            pidx = np.where(wwhere(fcad, 'plate*'))[0]
            for i, filt in enumerate(filters):
                hline = Span(location=0, dimension='width', line_color='black', line_alpha=.2)
                axs2[i].renderers.extend([hline])
                if len(pidx) > 0: axs2[i].scatter(mjds[pidx], fsn2[filt][pidx],  **plate_style)
                if len(bidx) > 0: axs2[i].scatter(mjds[bidx], fsn2[filt][bidx],  **bright_style)
                if len(didx) > 0: axs2[i].scatter(mjds[didx], fsn2[filt][didx],  **dark_style)
                if len(didx) > 0:
                    moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],fsn2[filt][didx]*900.0/exptime[didx])
                    line = axs2[i].line(moving_mjd, moving_avg, line_color="orange", line_width=1)
                    er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                    band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source = er, level='underlay', fill_alpha=.1, line_width=1, line_color='orange', fill_color='orange')
                    axs2[i].add_layout(band)
                if len(pidx) > 0:
                    moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],fsn2[filt][pidx]*900.0/exptime[pidx])
                    line = axs2[i].line(moving_mjd, moving_avg, line_color="blue", line_width=1)
                    er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                    band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source=er, level='underlay', fill_alpha=.1, line_width=1, line_color='blue', fill_color='blue')
                    axs2[i].add_layout(band)
            i = 3
            if len(pidx) > 0: axs2[3].scatter(mjds[pidx], see[pidx],  **plate_style)
            if len(bidx) > 0: axs2[3].scatter(mjds[bidx], see[bidx],  **bright_style)
            if len(didx) > 0: axs2[3].scatter(mjds[didx], see[didx],  **dark_style)
            if len(didx) > 0:
                moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],see[didx])
                line = axs2[i].line(moving_mjd, moving_avg, line_color="orange", line_width=1)
                er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source=er, level='underlay', fill_alpha=.1, line_width=1, line_color='orange', fill_color='orange')
                axs2[i].add_layout(band)
            if len(pidx) > 0:
                moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],see[pidx])
                line = axs2[i].line(moving_mjd, moving_avg, line_color="blue", line_width=1)
                er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source=er, level='underlay', fill_alpha=.1, line_width=1, line_color='blue', fill_color='blue')
                axs2[i].add_layout(band)

            if len(pidx) > 0:
                axs[3,0].scatter(mjds[pidx], see[pidx],  **plate_style)
                axs[3,1].scatter(mjds[pidx], see[pidx],  **plate_style)
                moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],see[pidx])
                line = axs[3,0].line(moving_mjd, moving_avg, line_color="blue", line_width=1)
                line = axs[3,1].line(moving_mjd, moving_avg, line_color="blue", line_width=1)
                er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source=er, level='underlay', fill_alpha=.1, line_width=1, line_color='blue', fill_color='blue')
                axs[3,0].add_layout(band)
                axs[3,1].add_layout(band)
            if len(bidx) > 0:
                axs[3,0].scatter(mjds[bidx], see[bidx], **bright_style)
                axs[3,1].scatter(mjds[bidx], see[bidx], **bright_style)
            if len(didx) > 0:
                axs[3,0].scatter(mjds[didx], see[didx],  **dark_style)
                axs[3,1].scatter(mjds[didx], see[didx],  **dark_style)
                moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],see[didx])
                line = axs[3,0].line(moving_mjd, moving_avg, line_color="orange", line_width=1)
                line = axs[3,1].line(moving_mjd, moving_avg, line_color="orange", line_width=1)
                er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source=er, level='underlay', fill_alpha=.1, line_width=1, line_color='orange', fill_color='orange')
                axs[3,0].add_layout(band)
                axs[3,1].add_layout(band)
        axs[3,0].legend.visible = True
        axs[3,0].legend.location = "top_left"
        axs[3,0].legend.label_text_font_size = '4pt'
        axs[3,0].legend.label_height = 3
        axs[3,0].legend.glyph_height = 3
        axs[3,0].legend.spacing = 0
        axs[3,0].legend.background_fill_alpha = .5
        axs[3,1].legend.visible = False

        if see is not None:
            ymin = np.nanmin(see)
            if ymin >= 0:
                axs2[i].y_range.start = -0.2
                axs[3,0].y_range.start = -0.2
                axs[3,1].y_range.start = -0.2
            ymax = np.nanmax(see)
            if ymax > 5.2:
                axs2[i].y_range.end = 5.2
                axs[3,0].y_range.end = 5.2
                axs[3,1].y_range.end = 5.2

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
                sos_data.to_csv(savename.replace('.html','_sos.csv').replace('BOSS_','.BOSS_'))
            else:
                _dir = ptt.dirname(savename)
                _name = ptt.basename(savename).replace('.html','_sos.csv').replace('BOSS_','.BOSS_')
                if ptt.exists(ptt.join(_dir, _name)):
                    sos_data = pd.read_csv(ptt.join(_dir, _name))
                else:
                    sos_data = None
                    if cron:splog.info(f'No OPSDB SOS data')
                    else: warnings.warn('No OPSDB SOS data')
            if cron: splog.info(f'Done Reading OPSDB Exposures')
        
        if sos_data is not None:
            for ccd in CCDs:
                sos_data[f'{ccd}_900'] = sos_data[ccd]/sos_data['exptime']*900
                axs3 = np.zeros([4],dtype=object)
                for i in range(len(axs3)):
                    ylim = None
                    xlab = 'MJD' if i == 1 else None
                    title = 'SOS SNR^2 Per Exposure (900s)' if i == 0 else None
                    if i < 2:
                        ylab = f'SOS SN2_{CCDs[i]}/(900s)'
                    elif i == 2:
                        ylab = 'SEEING50'
                    else:
                        ylab = 'Moon (%)'
                    
                    if i == 0:
                        try:
                            axs3[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab, y_range=ylim, title=title)
                        except:
                            axs3[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab,  title=title,sizing_mode='stretch_both')
                    else:
                        try:
                            axs3[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab, y_range=ylim, title=title, x_range = axs3[0].x_range)
                        except:
                            axs3[i] = figure(tools=tools,active_drag="box_zoom", x_axis_label=xlab,y_axis_label=ylab,  title=title,sizing_mode='stretch_both', x_range = axs3[0].x_range)
                    if i < 2:
                        axs3[i].add_tools(HoverTool(tooltips=[('MJD,SOS_SN2', '$x{int},$y{0.000}')]))
                    elif i == 2:
                        axs3[i].add_tools(HoverTool(tooltips=[('MJD,Seeing', '$x{int},$y{0.000}')]))
                    else:
                        axs3[i].add_tools(HoverTool(tooltips=[('MJD,Moon', '$x{int},$y{int}')]))
                    axs3[i].min_border_left = 75
            axs3, milestones = plot_milestone(obs, axs3, max_mjd, im=4, jm=1, html=True)
            dsidx = np.where(wwhere(sos_data['design_mode'].values, 'dark*'))[0]
            bsidx = np.where(wwhere(sos_data['design_mode'].values, 'bright*'))[0]
            for i, ccd in enumerate(CCDs):
                hline = Span(location=0, dimension='width', line_color='black', line_alpha=.2)
                axs3[i].renderers.extend([hline])
                if len(bsidx) > 0: axs3[i].scatter(sos_data['mjd'].values[bsidx], sos_data[f'{ccd}_900'].values[bsidx], **bright_style)
                if len(dsidx) > 0:
                    axs3[i].scatter(sos_data['mjd'].values[dsidx], sos_data[f'{ccd}_900'].values[dsidx], **dark_style)
                    moving_mjd, moving_avg, moving_16, moving_84 = runAvg(sos_data['mjd'].values[dsidx], sos_data[f'{ccd}_900'].values[dsidx])
                    line = axs3[i].line(moving_mjd, moving_avg, line_color="orange", line_width=1)
                    er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                    band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source=er, level='underlay', fill_alpha=.1, line_width=1, line_color='orange', fill_color='orange')
                    axs3[i].add_layout(band)
            i = 2
            if len(bidx) > 0: axs3[i].scatter(mjds[bidx], see[bidx],  **bright_style)
            if len(didx) > 0: axs3[i].scatter(mjds[didx], see[didx],  **dark_style)
            if len(didx) > 0:
                moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],see[didx])
                line = axs3[i].line(moving_mjd, moving_avg, line_color="orange", line_width=1)
                er = ColumnDataSource((pd.DataFrame({'moving_mjd':moving_mjd,'moving_16':moving_16,'moving_84':moving_84}).reset_index()))
                band = Band(base='moving_mjd', lower='moving_16', upper='moving_84', source=er, level='underlay', fill_alpha=.1, line_width=1, line_color='orange', fill_color='orange')
                axs3[i].add_layout(band)
            ymin = np.nanmin(see)
            if see is not None:
                if ymin >= 0:
                    axs3[i].y_range.start = -0.2
                ymax = np.nanmax(see)
                if ymax > 5.2:
                    axs3[i].y_range.end = 5.2
        
            i = 3
            axs3=axs3[0:3]
            axs3 = axs3.tolist()
            sosmjd = max(sos_data['mjd'])
        else:
            axs3= []
            sosmjd = ''
        
        l1 = layout(axs.tolist(), sizing_mode='stretch_both')
        if see is not None:
            l2 = layout(axs2.tolist(), sizing_mode='stretch_both')
        l3 = layout(axs3, sizing_mode='stretch_both')
        tab1 = Panel(child=l1,title='SpectroPhoto QA')
        if see is not None:
            tab2 = Panel(child=l2,title='SN2')
        tab3 = Panel(child=l3,title='SOS')
        if see is not None:
            tabs = Tabs(tabs=[tab1,tab2,tab3],sizing_mode='stretch_both')
        else:
            tabs = Tabs(tabs=[tab1,tab3], sizing_mode='stretch_both')
        title = Div(text="""<h2>SDSS-V """+obs+ """ BOSS QA: """+','.join(run2ds)+"""</h2>
                            <p>Spectro-photometric Plot Updated: """+all_mdate+ """<br>
                            Latest Full Pipeline MJD:"""+str(max(data['MJD']))+"""<br>
                            Latest SOS Pipeline MJD:"""+str(sosmjd) +"""<br>
                            Last Updated: """+time.ctime() +""" (MJD: """+str(int(Time.now().mjd))+""")</p>""",sizing_mode='stretch_width')
        save([title,tabs])
        print(savename)
        with open(savename) as f:
            soup = bs4.BeautifulSoup(f.read(),features="html5lib")
            style = soup.new_tag('style')
            soup.head.append(style)
            soup.select_one("style").append('.bk-logo {display:none !important;}')
            style=soup.new_tag('link rel="icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png"')
            soup.head.append(style)
        with open(savename, "w") as f:f.write(str(soup))

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


