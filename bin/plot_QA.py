#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')

import pandas as pd
from os import getenv
import numpy as np
import argparse
import os.path as ptt
from pydl.pydlutils.yanny import read_table_yanny
import time
from sdssdb.peewee.sdss5db.targetdb import database
import sdssdb
database.set_profile('pipelines')
from sdssdb.peewee.sdss5db.targetdb import Field, Cadence, DesignMode, Design, DesignToField
import os.path as ptt
from astropy.io import fits
from astropy.table import Table, unique, vstack
import re
from tqdm import tqdm
from glob import glob

try:
    daily_dir = getenv('DAILY_DIR')
except:
    daily_dir = ptt.join(getenv('HOME'),'daily')


def wwhere(array, value):
     value = value.replace('*','[\w]*')
     r = re.compile(value, re.IGNORECASE)
     ret = np.full(len(array), False)
     idx = [i for i, x in enumerate(array) if r.search(x)]
     ret[idx] = True
     return(ret)

def plot_QA(run2ds, test, mjds={}, obs='APO', testp='/test/sean/', clobber_lists=False, publish = False, epoch=False):
    fig, axs = plt.subplots(3,2, figsize=[12,6])

    if not ptt.exists(ptt.join(daily_dir,'etc', 'RM_fields')) or clobber_lists is True:
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
    else:
        with open(ptt.join(daily_dir,'etc', 'DarkMonitor_fields'), 'r') as f:
            monit_fields = np.array(f.readlines())
            monit_fields= np.char.replace(monit_fields, '\n','').astype(int).tolist()

    for i,run2d in enumerate(run2ds):    
        if test[i] is True: test_path = testp
        else: test_path = ''
        if not epoch:
            data=pd.read_csv(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,run2d,'spCalib_QA-'+run2d+'.csv'))
            mdate = time.ctime(ptt.getmtime(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,run2d,'spCalib_QA-'+run2d+'.csv')))
        else:
            data=pd.read_csv(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,run2d,'epoch','spCalib_QA-'+run2d+'.csv'))
            mdate = time.ctime(ptt.getmtime(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,run2d,'epoch','spCalib_QA-'+run2d+'.csv')))
        if mjds is not None:
            if run2d in mjds.keys():
                
                if mjds[run2d][0] is not None: data=data[data['MJD']>=mjds[run2d][0]]
                if mjds[run2d][1] is not None: data=data[data['MJD']<=mjds[run2d][1]]
#        print(data.to_string())
        if obs is not None:
            data=data[data['OBS'] == obs]
        max_mjd = np.max(data['MJD'].values)
        RM=data[data['FIELD'].isin(rm_fields)]
        Monit=data[data['FIELD'].isin(monit_fields)]
        
        milestones = read_table_yanny(ptt.join(daily_dir,'etc','fiber_milestones.par'), 'MILESTONE')
        milestones = milestones[np.where(milestones['obs'] == obs)]
        for row in milestones:
            for i in range(0,3):
                for j in range(0,2):
                    label = row['label'] if row['label'] != '' else None
                    if row['linestyle'] in ['-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted']:
                        ls = row['linestyle']
                    else:
                        ls = ( 0, tuple(np.array(tuple(row['linestyle'].replace(',','').replace(' ',''))).astype(int)))
                    if row['mjd'] > max_mjd + 1:
                        continue
                    axs[i,j].axvline(row['mjd'], color=row['color'], lw = .5, ls = ls, label=label)
        rm_style    = {'ls':'', 'marker':'.', 'color':'C1', 'alpha':.5, 'label':'RM Fields'}
        main_style  = {'ls':'', 'marker':'.', 'color':'C0', 'alpha':.2, 'label':'All Fields'}
        monit_style = {'ls':'', 'marker':'.', 'color':'C3', 'alpha':.5, 'label':'DarkMonitoring Fields'}
        ylim = [-0.10,0.10]
        ylim_r = [-0.001,0.15]

        filters = ['G','R','I']
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
    if not epoch:
        if len(run2ds) == 1:
            plt.savefig(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,run2d,'spCalib_QA-'+run2ds[0]+'-'+obs+'.png'), dpi=200)
        else: 
            plt.savefig(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,'spCalib_QA-'+'+'.join(run2ds)+'-'+obs+'.png'),dpi=200)
    else:
        if len(run2ds) == 1:
            plt.savefig(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,run2d,'epoch','spCalib_QA-'+run2ds[0]+'-'+obs+'.png'), dpi=200)
        else:
            plt.savefig(ptt.join(getenv("BOSS_SPECTRO_REDUX"),test_path,'spCalib_QA-epoch-'+'+'.join(run2ds)+'-'+obs+'.png'),dpi=200)
    plt.show()


    if len(run2ds) == 1:
        spall = None
        mjds = []
        mj = []
        fsn2g = []
        fsn2r = []
        fsn2i = []
        exptime = []
        fcad = []
        see = []

        ep = '' if not epoch else 'epoch'
        for i, row in tqdm(data.iterrows(), total=data.shape[0]):
            spallfile = glob(ptt.join(getenv("BOSS_SPECTRO_REDUX"), run2d, ep, 'spectra','full',str(row['FIELD']).zfill(6),str(row['MJD']),'spAll-'+str(row['FIELD']).zfill(6)+'-'+str(row['MJD'])+'.fits*'))
            if len(spallfile) == 0:
                continue
            else:
                spallfile = spallfile[0]
            try:
                spall = Table(fits.getdata(spallfile,1))
            except:
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

        fig, axs = plt.subplots(4,1, figsize=[12,8])
        didx = np.where(wwhere(fcad, 'dark*'))[0]
        bidx = np.where(wwhere(fcad, 'bright*'))[0]
        pidx = np.where(wwhere(fcad, 'plate*'))[0]
        for row in milestones:
            for i in range(0,4):
                label = row['label'] if row['label'] != '' else None
                if row['linestyle'] in ['-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted']:
                    ls = row['linestyle']
                else:
                    ls = ( 0, tuple(np.array(tuple(row['linestyle'].replace(',','').replace(' ',''))).astype(int)))
                if row['mjd'] > max(mj)+ 1:
                    continue
                axs[i].axvline(row['mjd'], color=row['color'], lw = .5, ls = ls, label=label)

        axs[0].grid(linestyle = ':', linewidth = 0.5)
        axs[0].minorticks_on() 
        axs[0].tick_params(axis='both', which='both', labelleft=True, labelright=True, labelbottom=False, labeltop=False, bottom=True, top=False, left=True, right=True)
        axs[0].tick_params(axis='x', which='both', bottom=True, top=True, direction='inout')
        axs[0].plot(mjds[pidx], fsn2g[pidx]*900.0/exptime[pidx], ls = '', marker= '.', ms = .5, color='C0', alpha=.5)
        axs[0].plot(mjds[bidx], fsn2g[bidx]*900.0/exptime[bidx], ls = '', marker= '.', ms = .5, color='C2', alpha=.5)
        axs[0].plot(mjds[didx], fsn2g[didx]*900.0/exptime[didx], ls = '', marker= '.', ms = .5, color='C1', alpha=.5)
        if len(didx) > 0:
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],fsn2g[didx]*900.0/exptime[didx])
            axs[0].plot(moving_mjd, moving_avg, color='C1', alpha=.5)
            axs[0].fill_between(moving_mjd,moving_16, moving_84, color='C1', alpha=.1)
        if len(pidx) > 0:
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],fsn2g[pidx]*900.0/exptime[pidx])
            axs[0].plot(moving_mjd, moving_avg, color='C0', alpha=.5)
            axs[0].fill_between(moving_mjd,moving_16, moving_84, color='C0', alpha=.1)
        axs[0].plot(np.nan, np.nan,label='plate',  ls = '', marker= '.', color='C0', ms = 2)
        axs[0].plot(np.nan, np.nan,label='bright', ls = '', marker= '.', color='C2', ms = 2)
        axs[0].plot(np.nan, np.nan,label='dark',   ls = '', marker= '.', color='C1', ms = 2)
        axs[0].legend(frameon=True,fontsize=3,loc=2)
        axs[0].set_xlabel('MJD')
        axs[0].set_ylabel('FieldSN2_G/(900s)')

        axs[1].grid(linestyle = ':', linewidth = 0.5)
        axs[1].minorticks_on() 
        axs[1].tick_params(axis='both', which='both', labelleft=True, labelright=True, labelbottom=False, labeltop=False, bottom=True, top=False, left=True, right=True)
        axs[1].tick_params(axis='x', which='both', bottom=True, top=True, direction='inout')
        axs[1].grid(linestyle = ':', linewidth = 0.5)
        axs[1].minorticks_on() 
        axs[1].plot(mjds[pidx], fsn2r[pidx]*900.0/exptime[pidx], ls = '', marker= '.', ms = .5,label='plate',  color='C0', alpha=.5)
        axs[1].plot(mjds[bidx], fsn2r[bidx]*900.0/exptime[bidx], ls = '', marker= '.', ms = .5,label='bright', color='C2', alpha=.5)
        axs[1].plot(mjds[didx], fsn2r[didx]*900.0/exptime[didx], ls = '', marker= '.', ms = .5,label='dark',   color='C1', alpha=.5)
        if len(didx) > 0:
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],fsn2r[didx]*900.0/exptime[didx])
            axs[1].plot(moving_mjd, moving_avg, color='C1', alpha=.5)
            axs[1].fill_between(moving_mjd,moving_16, moving_84, color='C1', alpha=.1)
        if len(pidx) > 0:
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],fsn2r[pidx]*900.0/exptime[pidx])
            axs[1].plot(moving_mjd, moving_avg, color='C0', alpha=.5)
            axs[1].fill_between(moving_mjd,moving_16, moving_84, color='C0', alpha=.1)
        axs[1].set_ylim(-.5,22)
        axs[1].set_xlabel('MJD')
        axs[1].set_ylabel('FieldSN2_R/(900s)')

        axs[2].grid(linestyle = ':', linewidth = 0.5)
        axs[2].minorticks_on() 
        axs[2].tick_params(axis='both', which='both', labelleft=True, labelright=True, labelbottom=True, labeltop=False, bottom=True, top=False, left=True, right=True)
        axs[2].tick_params(axis='x', which='both', bottom=True, top=True, direction='inout')
        axs[2].plot(mjds[pidx], fsn2i[pidx]*900.0/exptime[pidx], ls = '', marker= '.', ms = .5,label='plate',  color='C0', alpha=.5)
        axs[2].plot(mjds[bidx], fsn2i[bidx]*900.0/exptime[bidx], ls = '', marker= '.', ms = .5,label='bright', color='C2', alpha=.5)
        axs[2].plot(mjds[didx], fsn2i[didx]*900.0/exptime[didx], ls = '', marker= '.', ms = .5,label='dark',   color='C1', alpha=.5)
        if len(didx) > 0:
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],fsn2i[didx]*900.0/exptime[didx])
            axs[2].plot(moving_mjd, moving_avg, color='C1', alpha=.5)
            axs[2].fill_between(moving_mjd,moving_16, moving_84, color='C1', alpha=.1)
        if len(pidx) > 0:        
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],fsn2g[pidx]*900.0/exptime[pidx])
            axs[2].plot(moving_mjd, moving_avg, color='C0', alpha=.5)
            axs[2].fill_between(moving_mjd,moving_16, moving_84, color='C0', alpha=.1)
        axs[2].set_xlabel('MJD')
        axs[2].set_ylabel('FieldSN2_I/(900s)')
        axs[2].set_ylim(-.5,22)

        axs[3].grid(linestyle = ':', linewidth = 0.5)
        axs[3].minorticks_on()
        axs[3].tick_params(axis='both', which='both', labelleft=True, labelright=True, labelbottom=True, labeltop=False, bottom=True, top=False, left=True, right=True)
        axs[3].tick_params(axis='x', which='both', bottom=True, top=True, direction='inout')
        axs[3].plot(mjds[pidx], see[pidx], ls = '', marker= '.', ms = .5,label='plate',  color='C0', alpha=.5)
        axs[3].plot(mjds[bidx], see[bidx], ls = '', marker= '.', ms = .5,label='bright', color='C2', alpha=.5)
        axs[3].plot(mjds[didx], see[didx], ls = '', marker= '.', ms = .5,label='dark',   color='C1', alpha=.5)
        if len(didx) > 0:
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[didx],see[didx])
            axs[3].plot(moving_mjd, moving_avg, color='C1', alpha=.5)
            axs[3].fill_between(moving_mjd,moving_16, moving_84, color='C1', alpha=.1)
        if len(pidx) > 0:
            moving_mjd, moving_avg, moving_16, moving_84 = runAvg(mjds[pidx],see[pidx])
            axs[3].plot(moving_mjd, moving_avg, color='C0', alpha=.5)
            axs[3].fill_between(moving_mjd,moving_16, moving_84, color='C0', alpha=.1)
        axs[3].set_xlabel('MJD')
        axs[3].set_ylabel('SEEING50')
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
        if not epoch:
            plt.savefig(ptt.join(getenv("BOSS_SPECTRO_REDUX"), test_path, run2d, 'SN2-'+run2ds[0]+'-'+obs+'.png'), dpi=200)
        else:
            plt.savefig(ptt.join(getenv("BOSS_SPECTRO_REDUX"), test_path, run2d, 'epoch','SN2-'+run2ds[0]+'-'+obs+'.png'), dpi=200)
        plt.show()

def runAvg(mjd, val, ws=7):
    i = int(np.nanmin(mjd))
    moving_mjd = []
    moving_avg = []
    moving_16  = []
    moving_84  = []
    maxi = int(np.nanmax(mjd))#-ws+1
    while i < maxi:
        idx = np.where((mjd>=i) & (mjd < i+ws))[0]
        i=i+1
        if len(idx) == 0:
            moving_mjd.append(i+.5*ws)
            moving_avg.append(np.NaN)
            moving_16.append(np.NaN)
            moving_84.append(np.NaN)
        else:
            moving_mjd.append(np.nanmean(mjd[idx]))
            moving_avg.append(np.nanmean(val[idx]))
            moving_16.append(np.percentile(val[idx],16))
            moving_84.append(np.percentile(val[idx],84))
    return(moving_mjd, moving_avg, moving_16, moving_84)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot QA')
    parser.add_argument('-r','--run2d',    required=False, default = [getenv('RUN2D')], nargs="*", help='List of run2ds')
    parser.add_argument('-t','--test',     required=False, default = [False],    nargs="*", help='List of True/False test run2d (corresponding to run2d)')
    parser.add_argument('--test_path',     required=False, default='/test/sean/',help='test Run2d path modification')
    parser.add_argument('--mjds_low',      required=False, default = [None],     nargs="*", help='List of mjd lower limits - use None for no limit (corresponding to run2d)')
    parser.add_argument('--mjds_high',     required=False, default = [None],     nargs="*", help='List of mjd upper limits - use None for no limit (corresponding to run2d)')
    parser.add_argument('--clobber_lists', required=False, action='store_true',  help='Clobber list of fieldIDs')
    parser.add_argument('--lco',           required=False, action='store_true',  help='Flag for LCO vs APO')
    parser.add_argument('--publish',       required=False, action='store_true',  help='create publication version of plot')
    parser.add_argument('-e','--epoch',    required=False, action='store_true',  help='produce plots for epoch coadds')
    args = parser.parse_args()
    mjds={}
    
    
    for i, run2d in enumerate(args.run2d):
        if args.mjds_high[i] == 'None': args.mjds_high[i] = None 
        elif args.mjds_high[i] is None: continue
        else: args.mjds_high[i] = int(args.mjds_high[i])
        if args.mjds_low[i] == 'None': args.mjds_low[i] = None 
        elif args.mjds_low[i] is None: continue
        else: args.mjds_low[i] = int(args.mjds_low[i])
        mjds[run2d] = [args.mjds_low[i], args.mjds_high[i]]
    if args.lco is True: obs='LCO' 
    else: obs='APO'
    
    if len(args.test) != len(args.run2d): args.test.extend( [False] * (len(args.run2d)-len(args.test)))
    plot_QA(args.run2d, args.test, mjds=mjds, obs=obs, testp=args.test_path, clobber_lists=args.clobber_lists, publish= args.publish, epoch=args.epoch)
