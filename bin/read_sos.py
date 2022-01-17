#!/usr/bin/env python3       

from astropy.table import Table
import pandas as pd
from astropy.io import fits
import numpy as np
import warnings
import os.path as ptt
import matplotlib.pyplot as plt
import argparse
import sys
import glob
from os import mkdir,rename#,symlink
from shutil import copy

from pydl.pydlutils.trace import traceset2xy, TraceSet
from datetime import datetime
from pytz import timezone



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

def buildHTML(mjd, sos_dir='/data/boss/sos/'):
    figs = sorted(glob.glob(ptt.join(sos_dir,str(mjd).zfill(5),'summary_*')), key=ptt.getmtime, reverse=True)
    with open(ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'),'w') as f:
        f.write('<HTML>')
        f.write('<HEAD><TITLE>SOS plots for MJD '+str(mjd)+'</TITLE></HEAD>')
        reload_cmd="timerID=setTimeout('location.reload(true)',60000)"
        f.write('<BODY ONLOAD="'+reload_cmd+'">')
        f.write('<H2>SOS plots for MJD '+str(mjd)+' (<A HREF=../'+str(mjd)+'/logfile-'+str(mjd)+'.html>SOS Logs</A>)</H2>')
        # datetime object containing current date and time
        now = datetime.now(timezone('UTC'))
        f.write('<BR>This page last updated <B>'+now.strftime("%a %b %d %H:%M:%S %Y %Z")+'</B>.<P>')

        f.write('<TABLE BORDER=2>')
        plt_exps=[]
        for fig in figs:
            filenam=ptt.basename(fig)
            exp=str(int(fig.split('-')[1]))
            if exp in plt_exps: continue
            plt_exps.append(exp)
            ccd=fig.split('-')[-1].split('.')[0]
            b_fn ='../'+str(mjd).zfill(5)+'/summary_'+str(mjd).zfill(5)+'-'+exp.zfill(8)+'-b1.jpg'
            r_fn ='../'+str(mjd).zfill(5)+'/summary_'+str(mjd).zfill(5)+'-'+exp.zfill(8)+'-r1.jpg'
            b_fn_tn ='../'+str(mjd).zfill(5)+'/summary_'+str(mjd).zfill(5)+'-'+exp.zfill(8)+'-b1_wide.jpg'
            r_fn_tn ='../'+str(mjd).zfill(5)+'/summary_'+str(mjd).zfill(5)+'-'+exp.zfill(8)+'-r1_wide.jpg'          
            
            f.write('<TR><TD>'+exp+'<TD><A HREF='+b_fn+'> <IMG SRC='+b_fn_tn+' WIDTH=1200></A><br><A HREF='+r_fn+'> <IMG SRC='+r_fn_tn+' WIDTH=1200></A></TD>')
        f.write(' </TABLE></BODY></HTML>')
    copy(ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'),ptt.join(sos_dir,'combined','Summary_Current.tmp'))
    #symlink(ptt.join(sos_dir,str(mjd).zfill(5),'Summary_'+str(mjd).zfill(5)+'.html'), ptt.join(sos_dir,'combined','Summary_Current.tmp'))
    rename(ptt.join(sos_dir,'combined','Summary_Current.tmp'), ptt.join(sos_dir,'combined','Summary_Current.html'))

def Exp_summ(mjd, exposure, camera, sos_dir='/data/boss/sos/'):
    mjd=str(mjd)
    with fits.open(ptt.join(sos_dir,mjd,'logfile-'+mjd+'.fits')) as hdul:
        exp_log=Table(hdul[4].data)
    exp_log= exp_log[np.where((exp_log['EXPNUM']==exposure) & (exp_log['CAMERA']==camera))]

    PLUGFILE=exp_log['PLUGFILE'].value[0]
    CONFIGs=exp_log['CONFIG'].value[0]
    FIELD=exp_log['FIELD'].value[0]
    EXPNUM=exp_log['EXPNUM'].value[0]
    CAMERA=exp_log['CAMERA'].value[0]
    EXPTIME=exp_log['EXPTIME'].value[0]
    RAW_FLUX=exp_log['RAWFLUX'].value[0]
    RAW_FLUX_IVAR=exp_log['RAWFLUX_IVAR'].value[0]
    if 'DESIGNID' in exp_log.colnames: DESIGNID=exp_log['DESIGNID'].value[0]
    else: DESIGNID=''
    SN2=exp_log['SN2VECTOR'].value[0]
    mjd_exp=exp_log['TAI'].value[0]/(24.0*3600.0)

    wsetfile = glob.glob(ptt.join(sos_dir,str(mjd).zfill(5),'wset-'+str(mjd).zfill(5)+'-'+str(FIELD)+'-*-'+camera+'.fits'))[0]
    wset = fits.getdata(wsetfile, 1)
    xx, loglam = traceset2xy( TraceSet(wset))
    
    
    data = fits.getdata(ptt.join(sos_dir,mjd,'sci-'+str(CONFIGs).zfill(6)+'-'+camera+'-'+str(exposure).zfill(8)+'.fits'), ext = 0)
    if camera=='b1': 
        data=data[:,700:3500]
        wave = np.power(10.0,loglam[0])[700:3500]
    else: 
        data=data[:,300:3700]
        wave = np.power(10.0,loglam[0])[300:3700]
    exp_out=pd.DataFrame()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try: 
            with fits.open(ptt.join(sos_dir,mjd,'fibermap-'+str(FIELD)+'-'+camera+'.fits')) as hdul:
                extname='confSummary-'+str(CONFIGs)+'.par'
                head=hdul[extname].header
                plugmap=read_table(hdul[extname].data)

                if (bool(int(head['IS_DITHR'].split()[-1]))) or (head['PARENT_C'].split()[-1] != '-999'):
                    dithered=True
                    config_parent=head['PARENT_C'].split()[-1]
                    parent_extname='confSummary-'+str(config_parent)+'.par'
                    head_parent=hdul[0].header
                    plugmap_parent=read_table(hdul[1].data)
                    BOSSid = [k for k, i in enumerate(plugmap.FIBERTYPE.values) if 'BOSS' in i]

                    plugmap_parent=plugmap_parent.iloc[BOSSid]
                    plugmap_parent=plugmap_parent.sort_values('FIBERID', axis=0)
                else: dithered=False

            BOSSid = [k for k, i in enumerate(plugmap.FIBERTYPE.values) if 'BOSS' in i]

            plugmap=plugmap.iloc[BOSSid]#[plugmap.FIBERTYPE=='BOSS     ']
            plugmap=plugmap.sort_values('FIBERID', axis=0)

        except:
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
    exp_out["objtype"]=plugmap.OBJTYPE.values
    exp_out["ASSIGNED"]=plugmap.ASSIGNED.values
    exp_out["flux_g"]=plugmap.CALIBFLUX_1.values
    exp_out["flux_r"]=plugmap.CALIBFLUX_2.values
    exp_out["flux_i"]=plugmap.CALIBFLUX_3.values
    exp_out["flux_z"]=plugmap.CALIBFLUX_4.values

    exp_out["MAG_g"]=plugmap.MAG_1.values# +2.5*np.log10(2.085)
    exp_out["MAG_r"]=plugmap.MAG_2.values# +2.5*np.log10(2.085)
    exp_out["MAG_i"]=plugmap.MAG_3.values# +2.5*np.log10(2.116)
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

    exp_out_tab=Table.from_pandas(exp_out[['expid','exptime','mjd_obs','targetid','camera','target_ra','target_dec','fiber','objtype','flux_g','flux_r','flux_i','flux_z',
                               'hmag','spectroflux','spectroflux_ivar','delta_x_arcsec','delta_y_arcsec','xfocal','yfocal']])

    dither_path=ptt.join(sos_dir,mjd,'dither')
    if not ptt.isdir(dither_path): mkdir(dither_path)
    hdu=fits.table_to_hdu(exp_out_tab)
    hdu0 = fits.PrimaryHDU()
    hdu0.header['CONFIGID']=(int(CONFIGs), 'Configuration ID')
    hdulist = fits.HDUList([hdu0,hdu])
    hdulist.writeto(ptt.join(dither_path,'ditherBOSS-'+str(exposure).zfill(8)+'-'+camera+'-'+str(FIELD)+'.fits'),overwrite=True)
#    exp_out_tab.write(ptt.join(dither_path,'ditherBOSS-'+str(exposure).zfill(8)+'-'+camera+'-'+str(FIELD)+'.fits'),overwrite=True)
    return exp_out, wave, data, CONFIGs

def plot_exp(exp_out, wave, data, config, mjd, exp, ccd,log=True, sos_dir='/data/boss/sos/', wide=True, ref_data=None):

    def func1(t,a,b):
        t=np.power(10, 0.4*(22.5-t))
        return(np.power(a*t/np.sqrt(t+b),2))

    
    if wide is True:
        fig, axs = plt.subplot_mosaic([[0,1,2,3,4, 'z', 'z', 'z',]],
                                      constrained_layout=False,figsize=(21,3))
    else: 
        fig, axs = plt.subplot_mosaic([[0,1,2,3,4], ['z', 'z', 'z', 'z', 'z'], ['z', 'z', 'z', 'z', 'z']],
                                      constrained_layout=False,figsize=(12, 8))
    if ccd=='b1': filt='g'
    else: filt='i'
    
    skyid = [k for k, i in enumerate(exp_out.objtype.values) if 'SKY' in i]
    targid = [k for k, i in enumerate(exp_out.objtype.values) if 'SKY' not in i]
    skys=exp_out.iloc[skyid]
    targs=exp_out.iloc[targid]

    targs=targs[targs.MAG_g > -99]
    UnAssigned=targs[targs.ASSIGNED!=1]
    model_mags=np.arange(np.nanmin(targs['MAG_'+filt].values),np.nanmax(targs['MAG_'+filt].values), 0.001)
    model_mags=np.arange(np.nanmin(targs['MAG_'+filt].values),21.5, 0.001)

    #------------------------------------------
    axs[0].plot(exp_out.fiber, exp_out.spectroflux, ls='', marker='.',color='C0', mfc='none',label='obj')
    axs[0].plot(skys.fiber, skys.spectroflux, ls='', marker='.',color='C1', mfc='none',label='sky')
    if UnAssigned.shape[0]>0: axs[0].plot(UnAssigned.fiber, UnAssigned.spectroflux, ls='', marker='.', color='C2', mfc='none', label='UnAssigned')
    axs[0].legend(framealpha=.8)
    axs[0].set_xlabel('Fiberid')
    axs[0].set_ylabel('Flux')
    if log is True: axs[0].set_yscale('log')

    #------------------------------------------
    scale=np.nanmean((targs.exptime.values/900.0))
    axs[1].plot(targs['MAG_'+filt], targs.spectroflux, ls='', marker='.',color='C0', mfc='none')
    if ref_data is not None:
        axs[1].plot(ref_data['mag'],scale*ref_data['flux'], ls='', marker='.',color='C3', mfc='none')
    axs[1].set_xlabel('MAG_'+filt)
    axs[1].set_ylabel('Flux')
    scale=np.nanmean(targs.exptime.values/900.0)
    if ccd=='b1': 
        axs[1].plot(model_mags, scale*func1(model_mags,3.5704754364870914, 0.6684402741815383), 'k--',alpha=1)
        exp_out['dflux']=scale*func1(targs['MAG_'+filt],3.5704754364870914, 0.6684402741815383)-targs.spectroflux
    else: 
        axs[1].plot(model_mags, scale*func1(model_mags,3.826873867671731, -0.37855134101569876), 'k--',alpha=1)
        exp_out['dflux']=scale*func1(targs['MAG_'+filt],3.5704754364870914, 0.6684402741815383)-targs.spectroflux
    if log is True: axs[1].set_yscale('log')
    
    #------------------------------------------
    with np.errstate(invalid='ignore'):
        sc=axs[2].scatter(exp_out[exp_out.yfocal!=-999].xfocal, exp_out[exp_out.yfocal!=-999].yfocal, \
                          marker='.',c=np.log10(exp_out[exp_out.yfocal!=-999].dflux),cmap=plt.cm.jet,s=5,alpha=.8)
    plt.colorbar(sc,ax=axs[2],orientation='horizontal',label='log(ref_plate_flux-flux)',location='top')#,pad=0.01)
    axs[2].set_xlabel("x_focal")
    axs[2].set_ylabel("y_focal")
    axs[2].set_aspect('equal')
    
    
    #------------------------------------------
    axs[3].plot(targs.fiber, targs.SN2, ls='', marker='.',color='C0', mfc='none')
    axs[3].plot(skys.fiber, skys.SN2, ls='', marker='.',color='C1', mfc='none')
    if UnAssigned.shape[0]>0:axs[3].plot(UnAssigned.fiber, UnAssigned.SN2, ls='', marker='.', color='C2', mfc='none', label='UnAssigned')
    axs[3].set_xlabel('FIBERID')
    axs[3].set_ylabel('SN^2')
    if log is True: axs[3].set_yscale('log')

    #------------------------------------------
    scale=np.nanmean((targs.exptime.values/900.0))
    axs[4].plot(targs['MAG_'+filt], (targs.SN2), ls='', marker='.',color='C0', mfc='none')
    if ref_data is not None:
        axs[4].plot(ref_data['mag'],scale*ref_data['SNR'], ls='', marker='.',color='C3', mfc='none')
    axs[4].set_xlabel('MAG_'+filt)
    axs[4].set_ylabel('SN^2')
    if ccd=='b1':
        axs[4].plot(model_mags, scale*func1(model_mags, 3.673209821884937, 10.767534227684838), 'k--',alpha=1)
        exp_out['fsnr']=np.sqrt(exp_out.SN2)/np.sqrt(scale*func1(exp_out['MAG_'+filt],3.673209821884937, 10.767534227684838))
    else:
        axs[4].plot(model_mags, scale*func1(model_mags, 4.001601168006174, 26.750379730711874), 'k--',alpha=1)
        exp_out['fsnr']=np.sqrt(exp_out.SN2)/np.sqrt(scale*func1(exp_out['MAG_'+filt],4.001601168006174, 26.750379730711874))
    if log is True: axs[4].set_yscale('log')
    
    #------------------------------------------
    im=axs['z'].imshow(data, cmap='jet',resample=False,filternorm=False,aspect='auto',interpolation='none')
    axs['z'].set_ylabel("fiber")
    axs['z'].set_xlabel("Wavelength (Ang)")
    if ccd == 'b1':  x_label_list = [3500,4000,4500,5000,5500,6000,6500]
    else: x_label_list = [6000,6500,7000,7500,8000,8500,9000,9500,10000]
    axs['z'].set_xticks(find_nearest_indx(wave, x_label_list))
    axs['z'].set_xticklabels(x_label_list)
    
    fig.colorbar(im, orientation="vertical",label='flux',pad=0.01)

    fig.suptitle('MJD='+str(mjd)+'   EXP='+str(exp)+'   CCD='+ccd+'   CONFIG='+str(config)+'   FIELDID='+str(exp_out.iloc[0].fieldid)+'   DESIGNID='+str(exp_out.iloc[0].DESIGNID))
    fig.tight_layout()
    fig.show()
    if wide is True: fig.savefig(ptt.join(sos_dir,str(mjd), 'summary_'+str(mjd)+'-'+str(exp).zfill(8)+'-'+ccd+'_wide.jpg'))
    else: fig.savefig(ptt.join(sos_dir,str(mjd), 'summary_'+str(mjd)+'-'+str(exp).zfill(8)+'-'+ccd+'.jpg'), dpi=150)



def read_SOS(directory, mjd, exp=None, no_wide=False, ref_data=None):

    if exp is not None:
        exp=ptt.splitext(exp)[0]
        ccd=exp.split('-')[2]
        expNum=int(exp.split('-')[3])
        exp_out,wave,data,config=Exp_summ(mjd, expNum, ccd, sos_dir=directory)
        if no_wide is False: plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=True, ref_data=ref_data)
        plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=False, ref_data=ref_data)
        buildHTML(mjd,sos_dir=directory)
    else: 
        exps = sorted(glob.glob(ptt.join(directory, mjd, 'sci*.fits')), key=ptt.getmtime, reverse=True)
        for exp in exps:
            exp=ptt.splitext(exp)[0]
            ccd=exp.split('-')[2]
            expNum=int(exp.split('-')[3])
            exp_out,wave,data,config=Exp_summ(mjd, expNum, ccd, sos_dir=directory)
            if no_wide is False: plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=True, ref_data=ref_data)
            plot_exp(exp_out,wave,data,config,mjd, expNum, ccd, sos_dir=directory,wide=False, ref_data=ref_data)
            plt.close('all')
        buildHTML(mjd,sos_dir=directory)

if __name__ == '__main__' :
    parser = argparse.ArgumentParser(
        prog=ptt.basename(sys.argv[0]),
        description='Create Fiber info Summary')
    parser.add_argument('directory', type=str, help='SOS Directory')
    parser.add_argument('mjd', type=str, help='mjd')
    parser.add_argument('--exp', '-e',  type=str, help='Exposure Name', default=None)
    args = parser.parse_args()

    read_SOS(args.directory, args.mjd, args.exp)

