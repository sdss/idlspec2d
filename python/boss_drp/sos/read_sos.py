#!/usr/bin/env python3
from boss_drp.utils import find_nearest_indx
from boss_drp import idlspec2d_dir, favicon
from boss_drp.utils.hash import create_hash


from pydl.pydlutils.trace import traceset2xy, TraceSet
from pydl.pydlutils import yanny
from astropy.table import Table
from astropy.io import fits
import pandas as pd
import numpy as np
import warnings
from os import mkdir, rename, getenv, remove
import os.path as ptt
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
from time import sleep
import sys
import glob
from shutil import copy
from jinja2 import Template
import traceback
mpl.use('Agg')
plt.ioff()

def read_table(dataset):
    table=pd.DataFrame()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        for name in dataset.names:
            if len(dataset[name].shape)==1:
                data=dataset[name]
                if data.dtype.byteorder == '>':
                    data = data.byteswap().newbyteorder()
                table[name]=data
            else:
                for i in range(dataset[name].shape[1]):
                    data=dataset[name][:,i]
                    if data.dtype.byteorder == '>':
                        data = data.byteswap().newbyteorder()
                    table[name+'_'+str(i)]=data
    return(table)

def find_nearest_indx(array, value):
    indxs=np.zeros_like(value)
    array = np.asarray(array)
    for i, val in enumerate(value):
        indxs[i] = (np.abs(array - val)).argmin()
    return indxs

def find_confSummary(confid):
    SDSSCorepath = ptt.join(getenv('SDSSCORE_DIR'), getenv('OBSERVATORY').lower(), 'summary_files')
    SDSSCorepath = ptt.join(SDSSCorepath, str(int(np.floor(int(confid)/1000))).zfill(3)+'XXX')
    SDSSCorepath = ptt.join(SDSSCorepath, str(int(np.floor(int(confid)/100))).zfill(4)+'XX')
    if ptt.exists(ptt.join(SDSSCorepath,'confSummaryF-'+str(confid)+'.par')):
          confSummary = ptt.join(SDSSCorepath,'confSummaryF-'+str(confid)+'.par')
    elif ptt.exists(ptt.join(SDSSCorepath,'confSummary-'+str(confid)+'.par')):
          confSummary = ptt.join(SDSSCorepath,'confSummary-'+str(confid)+'.par')
    else: return(Table())
    confSummary=yanny.read_table_yanny(confSummary,'FIBERMAP')
    confSummary=confSummary[np.where(confSummary['fiberType'] == 'BOSS')]
    confSummary.sort('fiberId')
    return(confSummary)

def get_expid(filename):
    exp = ptt.basename(filename)
    exp = int(exp.split('-')[1])
    return(exp)

def buildHTML(mjd, sos_dir='/data/boss/sos/', nocopy=False, ccd=''):
    figs = sorted(glob.glob(ptt.join(sos_dir,str(mjd).zfill(5),'summary_*')), key=get_expid, reverse=True)
    
    template = ptt.join(idlspec2d_dir,'templates','html','SOS_Summary_template.html')
    out_file = ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'+ccd)

    try:
        now = datetime.now(datetime.UTC)
    except:
        now = datetime.utcnow()
    update = now.strftime("%a %b %d %H:%M:%S %Y %Z")
    
    
    if getenv('OBSERVATORY').lower() == 'apo': ccds= ['b1','r1']
    else: ccds= ['b2','r2']
    
    plt_exps=[]
    for fig in figs:
        filenam=ptt.basename(fig)
        exp=str(int(filenam.split('-')[1]))
        if exp.zfill(8) in plt_exps: continue
        plt_exps.append(exp.zfill(8))
    
    jinja2_opts = dict(MJD=mjd, date = now, UPDATE=update,
                        blue=ccds[0], red=ccds[1],
                        exps=plt_exps, favicon=favicon)
    
    with open(out_file, 'w', encoding="utf-8") as f:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            f.write(j2_template.render(jinja2_opts))

    try:
        rename(ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'+ccd), ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'))
    except:
        sleep(2)
        buildHTML(mjd, sos_dir=sos_dir, nocopy=nocopy, ccd = ccd)  

    if nocopy is False:
        try:
            copy(ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'),ptt.join(sos_dir,'combined','Summary_Current.tmp'))
            rename(ptt.join(sos_dir,'combined','Summary_Current.tmp'), ptt.join(sos_dir,'combined','Summary_Current.html'))
        except:
            sleep(2)
            try:
                copy(ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'),ptt.join(sos_dir,'combined','Summary_Current.tmp'))
                rename(ptt.join(sos_dir,'combined','Summary_Current.tmp'), ptt.join(sos_dir,'combined','Summary_Current.html'))
            except:
                print('Failed adding to combined directory')
                exit()

def Exp_summ(mjd, exposure, camera, sos_dir='/data/boss/sos/'):
    mjd=str(mjd)
    try: 
        with fits.open(ptt.join(sos_dir,mjd,'logfile-'+mjd+'.fits')) as hdul:
            exp_log=Table(hdul[4].data)
    except:
        sleep(2)
        #Try a secound time incase it was being written to
        try: 
            with fits.open(ptt.join(sos_dir,mjd,'logfile-'+mjd+'.fits')) as hdul:
                exp_log=Table(hdul[4].data)
        except:
            print('Failure opening '+'logfile-'+mjd+'.fits')
            exit()
    try:
        exp_log= exp_log[np.where((exp_log['EXPNUM']==exposure) & (exp_log['CAMERA']==camera))]
    except Exception as e:
        tb_str = traceback.format_exception(etype=type(e), value=e, tb=e.__traceback__)
        print("".join(tb_str))
        print('Invalid '+'logfile-'+mjd+'.fits'+' format')
        exit()
    if len(exp_log) == 0: exit()
    PLUGFILE=exp_log['PLUGFILE'].value[0]
    CONFIGs=exp_log['CONFIG'].value[0]
    FIELD=str(exp_log['FIELD'].value[0])
    FIELD=np.char.zfill(FIELD,6).tolist()
    EXPNUM=exp_log['EXPNUM'].value[0]
    CAMERA=exp_log['CAMERA'].value[0]
    EXPTIME=exp_log['EXPTIME'].value[0]
    RAW_FLUX=exp_log['RAWFLUX'].value[0]
    RAW_FLUX_IVAR=exp_log['RAWFLUX_IVAR'].value[0]
    if 'DESIGNID' in exp_log.colnames: DESIGNID=exp_log['DESIGNID'].value[0]
    else: DESIGNID=''
    SN2=exp_log['SN2VECTOR'].value[0]
    mjd_exp=exp_log['TAI'].value[0]/(24.0*3600.0)


    wsetfile = glob.glob(ptt.join(sos_dir,str(mjd).zfill(5),'wset-'+str(mjd).zfill(5)+'-'+str(FIELD)+'-*-'+camera+'.fits'))
    if len(wsetfile) == 0: FIELD=str(FIELD).zfill(6)
    wsetfile = glob.glob(ptt.join(sos_dir,str(mjd).zfill(5),'wset-'+str(mjd).zfill(5)+'-'+str(FIELD)+'-*-'+camera+'.fits'))
    if len(wsetfile) == 0:
        wsetfile = glob.glob(ptt.join(sos_dir,str(mjd).zfill(5),'wset-'+str(mjd).zfill(5)+'-*-*-'+camera+'.fits'))
        wsetfile.sort(key=ptt.getmtime)
        wsetfile.reverse()
    wset = fits.getdata(wsetfile[0], 1)
    xx, loglam = traceset2xy( TraceSet(wset))


    data = fits.getdata(ptt.join(sos_dir,mjd,'sci-'+str(CONFIGs).zfill(6)+'-'+camera+'-'+str(exposure).zfill(8)+'.fits'), ext = 0)
    hdr = fits.getheader(ptt.join(sos_dir,mjd,'sci-'+str(CONFIGs).zfill(6)+'-'+camera+'-'+str(exposure).zfill(8)+'.fits'), ext = 0)
    if (camera=='b1') or (camera=='b2'):
        data=data[:,700:3500]
        wave = np.power(10.0,loglam[0])[700:3500]
    else:
        data=data[:,300:3700]
        wave = np.power(10.0,loglam[0])[300:3700]
    exp_out=pd.DataFrame()
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      raw_mjd = mjd
      plugmap = None
      for mjd in [mjd, str(int(mjd)+1)]:
        if ptt.exists(ptt.join(sos_dir,mjd,'spfibermap-'+str(FIELD)+'-'+str(mjd)+'-'+camera+'.fits')):
            with fits.open(ptt.join(sos_dir,mjd,'spfibermap-'+str(FIELD)+'-'+str(mjd)+'-'+camera+'.fits')) as hdul:
                extname='confSummaryF-'+str(CONFIGs)+'.par'
                try: head=hdul[extname].header
                except: extname='confSummary-'+str(CONFIGs)+'.par'
                head=hdul[extname].header
                plugmap=read_table(hdul[extname].data)
                head = hdul['Summary'].data
                head = head[head['EXTNAME'] == extname][0]
                if (bool(int(head['IS_DITHERED'].split()[-1]))) or (head['PARENT_CONFIGURATION'].split()[-1] != '-999'):
                    dithered=True
                    yanny=False
                    config_parent=head['PARENT_CONFIGURATION'].split()[-1]
                    parent_extname='confSummaryF-'+str(config_parent)+'.par'
                    try: head_parent=hdul[parent_extname].header
                    except: 
                        parent_extname='confSummary-'+str(config_parent)+'.par'
                        try: 
                            head_parent = hdul['Summary'].data
                            head_parent = head_parent[head_parent['EXTNAME'] == parent_extname][0]
                        except:
                            yanny=True 
                            plugmap_parent=find_confSummary(config_parent)
                            #plugmap_parent=get_confSummary(config_parent, obs=None, sort=True, filter=True)
                    if yanny is False:
                        plugmap_parent=read_table(hdul[parent_extname].data)
                    else:
                        plugmap_parent.remove_column('mag')
                        plugmap_parent.rename_column('fiberType','FIBERTYPE')
                        plugmap_parent.rename_column('fiberId','FIBERID')
                        plugmap_parent.rename_column('ra','RA')
                        plugmap_parent.rename_column('dec','DEC')
                        plugmap_parent['FIBERTYPE'] = plugmap_parent['FIBERTYPE'].astype(str)
                        plugmap_parent = plugmap_parent.to_pandas()
                    BOSSid = [k for k, i in enumerate(plugmap_parent.FIBERTYPE.values) if 'BOSS' in i]

                    plugmap_parent=plugmap_parent.iloc[BOSSid]
                    plugmap_parent=plugmap_parent.sort_values('FIBERID', axis=0)
                else: dithered=False

            BOSSid = [k for k, i in enumerate(plugmap.FIBERTYPE.values) if 'BOSS' in i]

            plugmap=plugmap.iloc[BOSSid]#[plugmap.FIBERTYPE=='BOSS     ']
            plugmap=plugmap.sort_values('FIBERID', axis=0)
        elif ptt.exists(ptt.join(sos_dir,mjd,'fibermap-'+str(FIELD)+'-'+camera+'.fits')):
            with fits.open(ptt.join(sos_dir,mjd,'fibermap-'+str(FIELD)+'-'+camera+'.fits')) as hdul:
                extname='confSummaryF-'+str(CONFIGs)+'.par'
                try: head=hdul[extname].header
                except: extname='confSummary-'+str(CONFIGs)+'.par'
                head=hdul[extname].header
                plugmap=read_table(hdul[extname].data)

                if (bool(int(head['IS_DITHR'].split()[-1]))) or (head['PARENT_C'].split()[-1] != '-999'):
                    dithered=True
                    yanny=False
                    config_parent=head['PARENT_C'].split()[-1]
                    parent_extname='confSummaryF-'+str(config_parent)+'.par'
                    try: head_parent=hdul[parent_extname].header
                    except:
                        parent_extname='confSummary-'+str(config_parent)+'.par'
                        try: head_parent=hdul[parent_extname].header
                        except:
                            yanny=True
                            plugmap_parent=find_confSummary(config_parent)
                            #plugmap_parent=get_confSummary(config_parent, obs=None, sort=True, filter=True)
                    if yanny is False:
                        plugmap_parent=read_table(hdul[parent_extname].data)
                    else:
                        plugmap_parent.remove_column('mag')
                        plugmap_parent.rename_column('fiberType','FIBERTYPE')
                        plugmap_parent.rename_column('fiberId','FIBERID')
                        plugmap_parent.rename_column('ra','RA')
                        plugmap_parent.rename_column('dec','DEC')
                        plugmap_parent['FIBERTYPE'] = plugmap_parent['FIBERTYPE'].astype(str)
                        plugmap_parent = plugmap_parent.to_pandas()
                    BOSSid = [k for k, i in enumerate(plugmap_parent.FIBERTYPE.values) if 'BOSS' in i]

                    plugmap_parent=plugmap_parent.iloc[BOSSid]
                    plugmap_parent=plugmap_parent.sort_values('FIBERID', axis=0)
                else: dithered=False

            BOSSid = [k for k, i in enumerate(plugmap.FIBERTYPE.values) if 'BOSS' in i]

            plugmap=plugmap.iloc[BOSSid]#[plugmap.FIBERTYPE=='BOSS     ']
            plugmap=plugmap.sort_values('FIBERID', axis=0)

        elif ptt.exists(ptt.join(sos_dir,mjd,'fibermap-'+str(CONFIGs)+'-'+camera+'.fits')):
            with fits.open(ptt.join(sos_dir,mjd,'fibermap-'+str(CONFIGs)+'-'+camera+'.fits')) as hdul:
                head=hdul[0].header
                plugmap=read_table(hdul[1].data)

            BOSSid = [k for k, i in enumerate(plugmap.FIBERTYPE.values) if 'BOSS' in i]

            plugmap=plugmap.iloc[BOSSid]#[plugmap.FIBERTYPE=='BOSS     ']
            plugmap=plugmap.sort_values('FIBERID', axis=0)

            if (bool(int(head['is_dithe'].split()[-1]))) or (head['parent_c'].split()[-1] != '-999'):
                dithered=True
                config_parent=head['parent_c'].split()[-1]
                with fits.open(ptt.join(sos_dir,mjd,'fibermap-'+config_parent+'-'+camera+'.fits')) as hdul:
                    head_parent=hdul[0].header
                    plugmap_parent=read_table(hdul[1].data)
                BOSSid = [k for k, i in enumerate(plugmap.FIBERTYPE.values) if 'BOSS' in i]

                plugmap_parent=plugmap_parent.iloc[BOSSid]#[plugmap_parent.FIBERTYPE=='BOSS     ']
                plugmap_parent=plugmap_parent.sort_values('FIBERID', axis=0)

            else: dithered=False
        if plugmap is not None: break
    nfibs=len(plugmap.CATALOGID)
    if dithered:
        exp_out['dithered']=[dithered] * nfibs
        exp_out['parent']=[config_parent] * nfibs
    exp_out['CONFIG']=[CONFIGs] * nfibs
    exp_out['expid']=[EXPNUM] * nfibs
    exp_out['exptime']=[EXPTIME] * nfibs
    exp_out['fieldid']=[FIELD] * nfibs
    exp_out['DESIGNID']=[DESIGNID] * nfibs
    exp_out['mjd_obs']=[mjd_exp] * nfibs
    exp_out["targetid"]=plugmap.CATALOGID.values
    exp_out["camera"]=[camera] * nfibs
    if dithered:
        exp_out["target_ra"]=plugmap_parent.RA.values
        ra_parent=plugmap_parent.RA.values
        ra=plugmap.RA.values
        exp_out["target_dec"]=plugmap_parent.DEC.values
        dec_parent=plugmap_parent.DEC.values
        dec=plugmap.DEC.values
    else:
        exp_out["target_ra"]=plugmap.RA.values
        exp_out["target_dec"]=plugmap.DEC.values
    exp_out["fiber"]=plugmap.FIBERID.values
    try:
        exp_out['conffiberid']=plugmap.CONFFIBERID
    except:
        pass
    exp_out["objtype"]=plugmap.OBJTYPE.values
    exp_out["ASSIGNED"]=plugmap.ASSIGNED.values
    exp_out['valid']=plugmap.VALID.values
    exp_out['on_target']=plugmap.ON_TARGET.values
    exp_out["flux_g"]=plugmap.CALIBFLUX_1.values
    exp_out["flux_r"]=plugmap.CALIBFLUX_2.values
    exp_out["flux_i"]=plugmap.CALIBFLUX_3.values
    exp_out["flux_z"]=plugmap.CALIBFLUX_4.values

    exp_out["MAG_g"]=plugmap.MAG_1.values# +2.5*np.log10(2.085)
    exp_out["MAG_r"]=plugmap.MAG_2.values# +2.5*np.log10(2.085)
    exp_out["MAG_i"]=plugmap.MAG_3.values# +2.5*np.log10(2.116)

    if "CATDB_MAG_1" in plugmap:
        exp_out["CATDB_MAG_g"]=plugmap.CATDB_MAG_1.values
        exp_out["CATDB_MAG_r"]=plugmap.CATDB_MAG_2.values
        exp_out["CATDB_MAG_i"]=plugmap.CATDB_MAG_3.values
    else:  # deal with early commissioning data before implementation of optical_prov in readplugmap
        exp_out["CATDB_MAG_g"]=plugmap.MAG_1.values
        exp_out["CATDB_MAG_r"]=plugmap.MAG_2.values
        exp_out["CATDB_MAG_i"]=plugmap.MAG_3.values

    exp_out["SN2"]=SN2

    exp_out["hmag"]=plugmap.H_MAG.values
    exp_out["spectroflux"]=RAW_FLUX
    exp_out["spectroflux_ivar"]=RAW_FLUX_IVAR

    if dithered:
        exp_out['delta_x_arcsec'] = (ra-ra_parent)*np.cos(plugmap['DECFIELD'].values.astype(float)*np.pi/180.)*3600.
        exp_out['delta_y_arcsec'] = (dec-dec_parent)*3600.
    else:
        exp_out['delta_x_arcsec'] = 0.0
        exp_out['delta_y_arcsec'] = 0.0
    #exp_out["delta_x_arcsec"]=plugmap.DELTA_RA.values
    #exp_out["delta_y_arcsec"]=plugmap.DELTA_DEC.values
    exp_out["xfocal"]=plugmap.XFOCAL.values
    exp_out["yfocal"]=plugmap.YFOCAL.values
    
    for col in exp_out.columns:
        try:
            _data=exp_out[col].values
            if _data.dtype.byteorder == '>':
                _data = _data.byteswap().newbyteorder()
            exp_out[col]=_data
        except Exception as e:
            print(col, e)
            for i in range(exp_out[col].shape[1]):
                _data=exp_out[col][:,i].values
                if _data.dtype.byteorder == '>':
                    _data = _data.byteswap().newbyteorder()
                exp_out[name+'_'+str(i)]=[_data]
    
    exp_out_tab=Table.from_pandas(exp_out[['expid','exptime','mjd_obs','targetid','camera','target_ra','target_dec','fiber',
                                           'objtype','flux_g','flux_r','flux_i','flux_z',
                                           'MAG_g','MAG_r','MAG_i','CATDB_MAG_g','CATDB_MAG_r','CATDB_MAG_i','hmag',
                                           'spectroflux','spectroflux_ivar','delta_x_arcsec','delta_y_arcsec','xfocal','yfocal']])

    if 'conffiberid' in exp_out.columns:
        exp_out['conffiberid'] = exp_out['conffiberid'].astype(int)
        exp_out_tab=Table.from_pandas(exp_out[['expid','exptime','mjd_obs','targetid','camera','target_ra','target_dec','conffiberid',
                                               'objtype','flux_g','flux_r','flux_i','flux_z',
                                               'MAG_g','MAG_r','MAG_i','CATDB_MAG_g','CATDB_MAG_r','CATDB_MAG_i','hmag',
                                               'spectroflux','spectroflux_ivar','delta_x_arcsec','delta_y_arcsec','xfocal','yfocal']])
        exp_out_tab.rename_column('conffiberid','fiber')
    dither_path=ptt.join(sos_dir,mjd,'dither')
    try:
        if not ptt.isdir(dither_path): mkdir(dither_path)
    except:
        pass
    hdu=fits.table_to_hdu(exp_out_tab)
    hdu0 = fits.PrimaryHDU()
    hdu0.header['CONFIGID']=(int(CONFIGs), 'Configuration ID')
    hdu0.header['RA']=(hdr['RA'], 'RA of telescope boresight (deg)')
    try:hdu0.header['RADEG']=(hdr['RADEG'], 'RA of telescope pointing(deg)')
    except: pass
    hdu0.header['DEC']=(hdr['DEC'], 'DEC of telescope boresight (deg)')
    try:hdu0.header['DECDEG']=(hdr['DECDEG'], 'DEC of telescope pointing(deg)')
    except: pass
    try:hdu0.header['ROTPOS']=(hdr['ROTPOS'], 'Rotator request position (deg)')
    except: pass
    try:hdu0.header['AZ']=(hdr['AZ'], 'Azimuth axis pos. (approx, deg)')
    except: pass
    try: hdu0.header['ALT']=(hdr['ALT'], 'Altitude axis pos. (approx, deg)')
    except: pass
    try: hdu0.header['IPA']=(hdr['IPA'], 'Rotator axis pos. (approx, deg)')
    except: pass
    hdulist = fits.HDUList([hdu0,hdu])
    hdulist.writeto(ptt.join(dither_path,'ditherBOSS-'+str(exposure).zfill(8)+'-'+camera+'-'+str(FIELD)+'.fits'),overwrite=True)
    return exp_out, wave, data, CONFIGs

def plot_exp(exp_out, wave, data, config, mjd, exp, ccd,log=True, sos_dir='/data/boss/sos/', wide=True, ref_data=None):

    def func1(t,a,b):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            t=np.power(10, 0.4*(22.5-t))
            return(np.power(a*t/np.sqrt(t+b),2))


    if wide is True:
        fig, axs = plt.subplot_mosaic([[0,1,2,3,4, 'z', 'z', 'z',]],
                                      constrained_layout=False,figsize=(21,3))
    else:
        fig, axs = plt.subplot_mosaic([[0,1,2,3,4], ['z', 'z', 'z', 'z', 'z'], ['z', 'z', 'z', 'z', 'z']],
                                      constrained_layout=False,figsize=(12, 8))
    if (ccd=='b1') or (ccd=='b2'): filt='g'
    else: filt='i'

    skyid = [k for k, i in enumerate(exp_out.objtype.values) if 'SKY' in i]
    targid = [k for k, i in enumerate(exp_out.objtype.values) if 'SKY' not in i]
    skys=exp_out.iloc[skyid]
    targs=exp_out.iloc[targid]

    targs=targs[targs.MAG_g > -99]
    Assigned=targs[((targs.ASSIGNED==1) & (targs.valid==1))]
    model_mags=np.arange(np.nanmin(targs['MAG_'+filt].values),np.nanmax(targs['MAG_'+filt].values), 0.001)
    model_mags=np.arange(np.nanmin(targs['MAG_'+filt].values),21.5, 0.001)

    #------------------------------------------
    #axs[0].plot(exp_out.fiber, exp_out.spectroflux, ls='', marker='.',color='C0', mfc='none',label='obj')
    axs[0].plot(Assigned.fiber, Assigned.spectroflux, ls='', marker='.',color='C0', mfc='none',label='obj')
    axs[0].plot(skys.fiber, skys.spectroflux, ls='', marker='.',color='C1', mfc='none',label='sky')
    axs[0].legend(framealpha=.8)
    axs[0].set_xlabel('Fiberid')
    axs[0].set_ylabel('Flux')
    if log is True: axs[0].set_yscale('log')

    #------------------------------------------
    scale=np.nanmean((targs.exptime.values/900.0))
    axs[1].plot(Assigned['MAG_'+filt], Assigned.spectroflux, ls='', marker='.',color='C0', mfc='none')
    if ref_data is not None:
        axs[1].plot(ref_data['mag'],scale*ref_data['flux'], ls='', marker='.',color='C3', mfc='none')
    axs[1].set_xlabel('MAG_'+filt)
    axs[1].set_ylabel('Flux')
    scale=np.nanmean(Assigned.exptime.values/900.0)
    if ccd =='b1': fit = [3.5704754364870914, 0.6684402741815383]
    elif ccd== 'b2': fit = [3.32042298, 0.86944656]
    elif ccd== 'r1':fit = [3.826873867671731, -0.37855134101569876]
    elif ccd== 'r2':fit =[3.54210068,0.01387669]
    axs[1].plot(model_mags, scale*func1(model_mags,fit[0], fit[1]), 'k--',alpha=1)
    targs['dflux']=scale*func1(targs['MAG_'+filt],fit[0], fit[1])-targs.spectroflux
    targs['fflux']=scale*func1(targs['MAG_'+filt],fit[0], fit[1])/targs.spectroflux
    if log is True: axs[1].set_yscale('log')

    #------------------------------------------
    Assigned=targs[((targs.ASSIGNED==1) & (targs.valid==1))]
    with np.errstate(invalid='ignore'):
        sc=axs[2].scatter(targs[targs.yfocal!=-999].xfocal, targs[targs.yfocal!=-999].yfocal, \
                          marker='.',c=np.log10(targs[targs.yfocal!=-999].fflux),cmap=plt.cm.jet,s=5,alpha=.8)
    plt.colorbar(sc,ax=axs[2],orientation='horizontal',label='log(ref_plate_flux/flux)',location='top')#,pad=0.01)
    axs[2].set_xlabel("x_focal")
    axs[2].set_ylabel("y_focal")
    axs[2].set_aspect('equal')

    #------------------------------------------
    axs[3].plot(Assigned.fiber, Assigned.SN2, ls='', marker='.',color='C0', mfc='none')
    axs[3].plot(skys.fiber, skys.SN2, ls='', marker='.',color='C1', mfc='none')
    axs[3].set_xlabel('FIBERID')
    axs[3].set_ylabel('SN^2')
    if log is True: axs[3].set_yscale('log')

    #------------------------------------------
    scale=np.nanmean((targs.exptime.values/900.0))
    axs[4].plot(Assigned['MAG_'+filt], (Assigned.SN2), ls='', marker='.',color='C0', mfc='none')
    if ref_data is not None:
        axs[4].plot(ref_data['mag'],scale*ref_data['SNR'], ls='', marker='.',color='C3', mfc='none')
    axs[4].set_xlabel('MAG_'+filt)
    axs[4].set_ylabel('SN^2')
    if (ccd=='b1'): fit = [3.673209821884937, 10.767534227684838]
    elif (ccd=='b2'):fit = [3.31581073, 15.33823938]
    elif (ccd=='r1'):fit = [4.001601168006174, 26.750379730711874]
    elif (ccd=='r2'):fit = [3.43656542, 16.24230991]
    axs[4].plot(model_mags, scale*func1(model_mags, fit[0],fit[1]), 'k--',alpha=1)
    targs['fsnr']=np.sqrt(targs.SN2)/np.sqrt(scale*func1(targs['MAG_'+filt],fit[0], fit[1]))
    if log is True: axs[4].set_yscale('log')

    #------------------------------------------
    im=axs['z'].imshow(data, cmap='jet',resample=False,filternorm=False,aspect='auto',interpolation='none')
    axs['z'].set_ylabel("fiber")
    axs['z'].set_xlabel("Wavelength (Ang)")
    if (ccd=='b1') or (ccd=='b2'):  x_label_list = [3500,4000,4500,5000,5500,6000,6500]
    else: x_label_list = [6000,6500,7000,7500,8000,8500,9000,9500,10000]
    axs['z'].set_xticks(find_nearest_indx(wave, x_label_list))
    axs['z'].set_xticklabels(x_label_list)

    fig.colorbar(im, orientation="vertical",label='flux',pad=0.01)

    fig.suptitle('MJD='+str(mjd)+'   EXP='+str(exp)+'   CCD='+ccd+'   CONFIG='+str(config)+'   FIELDID='+str(exp_out.iloc[0].fieldid)+'   DESIGNID='+str(exp_out.iloc[0].DESIGNID))
    fig.tight_layout()
    #fig.show()
    if wide is True: fig.savefig(ptt.join(sos_dir,str(mjd), 'summary_'+str(mjd)+'-'+str(exp).zfill(8)+'-'+ccd+'_wide.jpg'))
    else: fig.savefig(ptt.join(sos_dir,str(mjd), 'summary_'+str(mjd)+'-'+str(exp).zfill(8)+'-'+ccd+'.jpg'), dpi=150)



def read_SOS(directory, mjd, exp=None, no_wide=False, ref_data=None, nocopy=False, update_hash=False):

    if exp is not None:
        exp=ptt.splitext(exp)[0]
        ccd=exp.split('-')[2]
        expNum=int(exp.split('-')[3])
        exp_out,wave,data,config=Exp_summ(mjd, expNum, ccd, sos_dir=directory)
        if no_wide is False: plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=True, ref_data=ref_data)
        plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=False, ref_data=ref_data)
        buildHTML(mjd,sos_dir=directory,nocopy=nocopy, ccd = ccd)
        
        html_file = ptt.join(directory,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'+ccd)
        if ptt.exists(html_file):
            remove(html_file)
    else:
        exps = sorted(glob.glob(ptt.join(directory, mjd, 'sci*.fits')), key=ptt.getmtime, reverse=False)
        for exp in exps:
            exp=ptt.splitext(exp)[0]
            ccd=exp.split('-')[2]
            expNum=int(exp.split('-')[3])
            exp_out,wave,data,config=Exp_summ(mjd, expNum, ccd, sos_dir=directory)
            if no_wide is False: plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=True, ref_data=ref_data)
            plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=False, ref_data=ref_data)
            plt.close('all')
        buildHTML(mjd,sos_dir=directory,nocopy=nocopy)



    if update_hash:
        test = create_hash(ptt.join(directory,mjd))
        if test:
            print("\nsha1sum is locked")
