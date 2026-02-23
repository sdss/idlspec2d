#!/usr/bin/env python3
from boss_drp import favicon, idlspec2d_dir
from boss_drp.prep.GetconfSummary import  get_confSummary
from boss_drp.utils import chpc2html
from boss_drp.Flatlib.plot import (plot_flat, plot_raw, get_raw, build_savename,
                                    plot_thruput_v_sextant, read_fiberAssignments,
                                    plot_thruput_timeseries)
from jinja2 import Template

from astropy.io import fits,ascii
from astropy.table import Table
from astropy.table import unique as tabunique
import argparse
from os import getenv, makedirs, remove
from glob import glob
import os.path as ptt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.colors as pc

import numpy as np
import pandas as pd
import random
from matplotlib import pyplot as plt
from datetime import date, datetime
from astropy.time import Time
from pydl.pydlutils import yanny
import warnings
import bs4
from tqdm import tqdm
today=datetime.now()


def ID_lowFiber(filename, save_dir, med, confFiberid, threshold=0.8, clobber=False):
    if np.isnan(med).all():
        return()
    makedirs(save_dir, exist_ok=True)
    flagged = np.where(med < threshold)
    if len(med.shape) == 1:
        flagged=flagged[0]
    else:
        flagged=flagged[1]
        flagged=np.unique(flagged)
    if len(confFiberid.shape)!=1:
        confFiberid = confFiberid[0]
        
    #flagged = flagged + 1
    flagged = np.column_stack((flagged + 1, confFiberid[flagged]))
    filename=ptt.basename(filename)
    filename=ptt.splitext(ptt.splitext(filename)[0])[0]
    filename=filename+'_lt'+str(threshold)+'.txt'
    if (not ptt.exists(ptt.join(save_dir,filename))):
        np.savetxt(ptt.join(save_dir,filename), flagged, fmt='%d', delimiter=',', header='Trace,Slit')
    elif clobber is True:
        np.savetxt(ptt.join(save_dir,filename), flagged, fmt='%d', delimiter=',', header='Trace,Slit')
    return



def find_confSummary(confid, obs):
    confSummary = get_confSummary(confid, obs=obs, no_remote=True, sort= True, filter=True)
    return(confSummary)

def create_att_path(obs,mjd, expsids,ccd=None, ftype='sdR', subdir=''):
    #fpath=np.char.add(obs,'/'+subdir)
    fpath=np.char.add(subdir,mjd)
    fpath=np.char.add(fpath,'/'+ftype+'-')
    if ccd is not None: fpath=np.char.add(fpath,ccd+'-')
    fpath=np.char.add(fpath, expsids)
    return(fpath)

def make_beta_plots(Allflat,directory, version, obsf, meds_b, meds_r, qbad_b, qbad_r, beta_b, beta_r, fiberids, beta, dbeta, subf):
    now =today.strftime("Last Updated: %a %b %d %H:%M:%S %Y (MJD: "+'{:.3f}'.format(Time(today).mjd)+")")
    outtitle = f'beta {beta}-{beta+dbeta}'
    output_dir=ptt.join(directory, obsf, 'beta',subf)
    outfile=ptt.join(output_dir,f'beta{str(beta).zfill(3)}-{str(beta+dbeta).zfill(3)}.html')
    title = '\n'.join([f"<h2>SDSS-V {obsf.upper()} BOSS Flat Libray Analysis:{version}</h2>",
                       f"<h2>{beta} < BETA < {beta+dbeta}</h2>",
                       f"<p>{now}</p>"])
    mask_b = np.full_like(meds_b, np.NaN)
    mask_b[np.where((beta_b >= beta) & (beta_b <=beta+dbeta))] = 1
    mask_b[:,np.where(np.nansum(mask_b,axis=0)<=1)[0]] = np.NaN

    mask_r = np.full_like(meds_r, np.NaN)
    mask_r[np.where((beta_r >= beta) & (beta_r <=beta+dbeta))] = 1
    mask_r[:,np.where(np.nansum(mask_r,axis=0)<=1)[0]] = np.NaN

    if not ((np.nansum(mask_b) > 0) or (np.nansum(mask_r) > 0)):
        if ptt.exists(outfile):
            remove(outfile)
        return
    makedirs(output_dir,exist_ok=True)

    med_plot(Allflat,directory,outfile, version, obsf, meds_b*mask_b, meds_r*mask_r,
             qbad_b, qbad_r, fiberids, title=title, name=outtitle,
             mask_b=mask_b, mask_r=mask_r)
             
             
def make_pos_plots(Allflat,directory, version, obs, meds_b, meds_r, qbad_b, qbad_r, fiberids, drot, dalt, alt, rot, subf):
    now =today.strftime("Last Updated: %a %b %d %H:%M:%S %Y (MJD: "+'{:.3f}'.format(Time(today).mjd)+")")
    if drot is None:
        output_dir=ptt.join(directory, obs, 'alt',subf)
        indx = np.where((Allflat['ALT'].value >= alt) & (Allflat['ALT'].value <= alt+dalt))[0]
        outfile=ptt.join(output_dir,f'alt{str(alt).zfill(3)}-{str(alt+dalt).zfill(3)}.html')
        outtitle=f'alt{alt}-{alt+dalt}'
        title = '\n'.join([f"<h2>SDSS-V {obs.upper()} BOSS Flat Libray Analysis:{version}</h2>",
                           f"<h2>{alt} < ALT < {alt+dalt}</h2>",
                           f"<p>{now}</p>"])
    elif dalt is None:
        output_dir=ptt.join(directory, obs, 'rot',subf)
        indx = np.where((Allflat['ROT'].value >= rot) & (Allflat['ROT'].value <= rot+drot))[0]
        outfile=ptt.join(output_dir,f'rot{str(rot).zfill(3)}-{str(rot+drot).zfill(3)}.html')
        outtitle=f'rot{rot}-{rot+drot}'
        title = '\n'.join([f"<h2>SDSS-V {obs.upper()} BOSS Flat Libray Analysis:{version}</h2>",
                           f"<h2>{rot} < ROT < {rot+drot}</h2>",
                           f"<p>{now}</p>"])
    else:
        output_dir=ptt.join(directory, obs, 'tele_pos',subf)
        indx = np.where((Allflat['ALT'].value >= alt) & (Allflat['ALT'].value <= alt+dalt) &
                        (Allflat['ROT'].value >= rot) & (Allflat['ROT'].value <= rot+drot))[0]
        outfile=ptt.join(output_dir,f'alt{str(alt).zfill(3)}-{str(alt+dalt).zfill(3)}_rot{str(rot).zfill(3)}-{str(rot+drot).zfill(3)}.html')
        outtitle=f'alt{alt}-{alt+dalt}_rot{rot}-{rot+drot}'
        title = '\n'.join([f"<h2>SDSS-V {obs.upper()} BOSS Flat Libray Analysis:{version}</h2>",
                           f"<h2>{alt} < ALT < {alt+dalt}; {rot} < ROT < {rot+drot}</h2>",
                           f"<p>{now}</p>"])
    if len(indx) == 0:
        if ptt.exists(outfile):
            remove(outfile)
        return
    makedirs(output_dir,exist_ok=True)
    med_plot(Allflat[indx],directory,outfile, version, obs, meds_b[indx], meds_r[indx],
             qbad_b[indx], qbad_r[indx], fiberids[indx], title=title, name=outtitle)

def position_plots(Allflat,directory, version, obs, meds_b, meds_r, qbad_b, qbad_r,beta_b, beta_r, fiberids, plates=False):
    makedirs(ptt.join(directory, obs, 'tele_pos'),exist_ok=True)
    makedirs(ptt.join(directory, obs, 'Rot'),     exist_ok=True)
    makedirs(ptt.join(directory, obs, 'alt'),     exist_ok=True)
    if plates is False: makedirs(ptt.join(directory, obs, 'beta'),exist_ok=True)
    fibers=fiberids#np.tile(np.arange(1,501),(len(qbad_b),1))

    tqdm.write('Making drot=dalt plots')
    for step in [5,10,15]:
        for alt in np.arange(25,90,step):
            for rot in np.arange(0, 360, step):
                make_pos_plots(Allflat,directory, version, obs, meds_b, meds_r, qbad_b, qbad_r, fibers,  step, step, alt, rot, str(step))

    tqdm.write('Making drot or dalt plots')
    for step in [5,10,15]:
        for deg in np.arange(25,90,step):
            make_pos_plots(Allflat, directory, version, obs, meds_b, meds_r, qbad_b, qbad_r, fibers, None, step, deg, None, str(step))
            make_pos_plots(Allflat, directory, version, obs, meds_b, meds_r, qbad_b, qbad_r, fibers, step, None, None, rot, str(step))
    tqdm.write('Making dBeta plots')
    if plates is False:
        for step in [5,10,15]:
            for beta in np.arange(0,180,step):
                make_beta_plots(Allflat,directory, version, obs, meds_b, meds_r, qbad_b, qbad_r, beta_b, beta_r, fibers, beta, step, str(step))

def hex_to_rgb_tuple(hex_color):
    # Convert hex to RGB tuple
    rgb_tuple = tuple(int(hex_color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
    return f'rgb({rgb_tuple[0]},{rgb_tuple[1]},{rgb_tuple[2]})'

def generate_discrete_sequence(base_scale, num_colors):
    # Convert all base scale colors to RGB strings
    base_scale_rgb = [hex_to_rgb_tuple(color) for color in base_scale]

    # Create a linear space for interpolating colors across the entire base scale
    t_values = np.linspace(0, 1, num_colors)
    
    # Generate the discrete color sequence
    discrete_sequence = []
    for t in t_values:
        # Find the appropriate position in the base scale
        scaled_t = t * (len(base_scale_rgb) - 1)
        idx = int(scaled_t)
        fractional_t = scaled_t - idx
        
        # Handle edge case where idx is at the last color
        if idx >= len(base_scale_rgb) - 1:
            discrete_sequence.append(base_scale_rgb[-1])
        else:
            color = pc.find_intermediate_color(base_scale_rgb[idx], base_scale_rgb[idx + 1], fractional_t, colortype='rgb')
            discrete_sequence.append(color)
    
    return discrete_sequence
    
def col2html(row, raw=False, ccd = 'Blue', directory=''):
    bsd = 'BOSS_SPECTRO_DATA_N' if row['Obs'] == 'apo' else 'BOSS_SPECTRO_DATA_S'
    bsd = chpc2html(getenv(bsd))
    if raw:
        template = ('<a href="./{value}.png" target="_blank">PNG</a> '+
                    '<a href="{bsd}/{value}.fit.gz" '+
                    'target="_blank">FITs</a>')
    else:
        template = ('<a href="./{value}.png" target="_blank">PNG</a> '+
                    '<a href="{directory}/{obs}/{value}.fits.gz" target="_blank">FITs</a> '+
                    '<a href="./{value}_trans.png" target="_blank">Trans</a>')
    col = f'{ccd} Raw' if raw else f'{ccd} Reduced'
    return(template.format(value=row[col], directory = chpc2html(directory),
                            bsd=bsd, obs= row['Obs']))

def med_plot(Allflat,directory,filename, ver, obsf, meds_b, meds_r, qbad_b, qbad_r,
            fiberids, title=None, name=None, mask_b=None, mask_r=None):
    ccds = ['b1','r1'] if obsf == 'apo' else ['b2','r2']
    umjs = np.unique(Allflat['MJD'])
    discrete_palette = generate_discrete_sequence(pc.sequential.Viridis, len(umjs))
    random.shuffle(discrete_palette)
    
    obs = Allflat['OBS'].value.tolist()
    mjd =Allflat['MJD'].value.astype(str).tolist()
    expsids = np.char.zfill(Allflat['EXP'].value.astype(str),8)
    red_reduced=create_att_path(obs,mjd,expsids,ccd=ccds[1],ftype='spFlat')
    red_raw=create_att_path(obs,mjd,expsids,ccd=ccds[1],ftype='sdR')
    blue_reduced=create_att_path(obs,mjd,expsids,ccd=ccds[0],ftype='spFlat')
    blue_raw=create_att_path(obs,mjd,expsids,ccd=ccds[0],ftype='sdR')
    red_trans=create_att_path(obs,mjd,expsids,ccd=ccds[1], ftype='spFlat', subdir='transmission/')
    blue_trans=create_att_path(obs,mjd,expsids,ccd=ccds[0], ftype='spFlat', subdir='transmission/')

    meta = pd.DataFrame({'exp':np.char.zfill(Allflat['EXP'].value.astype(str),8),
                         'MJD':Allflat['MJD'].value.astype(str).tolist(),
                         'TAI':Allflat['TAI'].value.tolist(),
                         'qbad_b':qbad_b,
                         'qbad_r':qbad_r,
                         'airmass':Allflat['AIRMASS'].value.tolist(),
                         'Alt':Allflat['ALT'].value.tolist(),
                         'rot':Allflat['ROT'].value.tolist(),
                         'Obs':Allflat['OBS'].value.tolist(),
                         'Blue Raw':blue_raw,
                         'Blue Reduced':blue_reduced,
                         'Red Raw':red_raw,
                         'Red Reduced':red_reduced})

    meta['Blue Raw'] = meta.apply(lambda row: col2html(row, raw = True, ccd = 'Blue', directory=directory), axis=1)
    meta['Blue Reduced'] = meta.apply(lambda row: col2html(row, ccd = 'Blue', directory=directory), axis=1)
    meta['Red Raw'] = meta.apply(lambda row: col2html(row, raw = True, ccd = 'Red', directory=directory), axis=1)
    meta['Red Reduced'] = meta.apply(lambda row: col2html(row, raw = False, ccd = 'Red', directory=directory), axis=1)
    
    hr = '<br>'.join(["exp: {exp}, qbad: ({qbad_b},{qbad_r})",
                      "MJD: {MJD}, TAI: {TAI}",
                      "Airmass: {airmass}, Alt: {Alt}, Rot: {rot}",
                      "Observatory: {Obs}"])
    ht = [hr.format(**row.to_dict()) for i,row in meta.iterrows()]
    figs = {}
    for type in ['all','select']:
        figs[type] = make_subplots(rows=2, cols=1, shared_xaxes=True, shared_yaxes=False,vertical_spacing=0.02)
        groups = []
        for i in range(len(qbad_b)):
            vis = True
            if type == 'select':
                vis = True if Allflat['MJD'][i] == max(umjs) else 'legendonly'
                    
            color = discrete_palette[np.where(umjs == Allflat['MJD'][i])[0][0]]
            if str(Allflat['MJD'][i]) not in groups:
                showlegend = True
                groups.append(str(Allflat['MJD'][i]))
            else: showlegend = False
            params = {'mode':'lines','name': str(Allflat['MJD'][i]),
                      'legendgroup':str(Allflat['MJD'][i]), 'visible':vis,
                      'line':dict(color=color,dash='solid', width=.5,)}
            params['text']=[ht[0]]*len(fiberids[i])
            params['hovertemplate']='<br>Fiber:%{x}, Throughput:%{y}<br>%{text}'
            figs[type].add_trace(go.Scatter(x= fiberids[i], y=meds_b[i],
                                            showlegend=showlegend,**params), row=1, col=1)
            figs[type].add_trace(go.Scatter(x= fiberids[i], y=meds_r[i], showlegend=False,
                                            **params), row=2, col=1)
        figs[type].update_layout(autosize=True, height=None, margin=dict(l=0,r=0,t=20,b=0))
        figs[type].update_xaxes(title='Trace Fiber', row=2, col=1)
        figs[type].update_yaxes(title='Med Flat Blue', row=1, col=1)
        figs[type].update_yaxes(title='Med Flat Red', row=2, col=1)
        
    if mask_b is None:
        mask_b=np.isfinite(meds_b)
    if mask_r is None:
        mask_r=np.isfinite(meds_r)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ys_b = {'div':(np.nanmax(meds_b,axis=0) - np.nanmin(meds_b,axis=0))/np.nanmean(meds_b,axis=0)*100,
                'std':np.nanstd(meds_b,ddof=1,axis=0),
                'mad':np.nanmedian(np.absolute(meds_b - np.nanmedian(meds_b, axis=0)),axis=0),
                'count':np.nansum(mask_b,axis=0)}
        ys_r = {'div':(np.nanmax(meds_r,axis=0) - np.nanmin(meds_r,axis=0))/np.nanmean(meds_r,axis=0)*100,
                'std':np.nanstd(meds_r,ddof=1,axis=0),
                'mad':np.nanmedian(np.absolute(meds_r - np.nanmedian(meds_r, axis=0)),axis=0),
                'count':np.nansum(mask_r,axis=0)}
    labels ={'div':'% Max Dev','std':'STD','mad':'Med Abs Deviation','count':'N Flats'}
    for type in ys_b.keys():
        figs[type] = make_subplots(rows=2, cols=1, shared_xaxes=True, shared_yaxes=False,vertical_spacing=0.02)

        params = {'mode':'lines','name': type, 'legendgroup':type,
                  'line':dict(color='black',dash='solid', width=.5)}
        params['text']=[ht[0]]*len(fiberids[i])
        params['hovertemplate']='<br>Fiber:%{x}, '+labels[type]+':%{y}<br>%{text}'

        figs[type].add_trace(go.Scatter(x= fiberids[0], y=ys_b[type], showlegend=False,**params), row=1, col=1)
        figs[type].add_trace(go.Scatter(x= fiberids[0], y=ys_r[type], showlegend=False,**params), row=2, col=1)
        figs[type].update_layout(autosize=True, height=None, margin=dict(l=0,r=0,t=20,b=0))
        figs[type].update_xaxes(title='Trace Fiber', row=2, col=1)
        figs[type].update_yaxes(title=f'{labels[type]} Blue', row=1, col=1)
        figs[type].update_yaxes(title=f'{labels[type]} Red', row=2, col=1)

    if title is None:
        now =today.strftime("Last Updated: %a %b %d %H:%M:%S %Y (MJD: "+'{:.3f}'.format(Time(today).mjd)+")")
        title = '\n'.join([f"<h2>SDSS-V {obsf.upper()} BOSS Flat Libray Analysis:{ver}</h2>",
                           f"<p>{now}</p>",
                           f"<p>N Flats Included: {str(np.nansum(mask_b,axis=0)[0])}</p>"])
 
    if name is None:
        name = f"{obsf.upper()} BOSS Flats: {ver}"
    
    fignames = {'all':f'BOSSFlats_{obsf.upper()}_{ver}',
                'select':f'BOSSFlats_{obsf.upper()}_{ver}_select',
                'div':f'BOSSFlats_{obsf.upper()}_{ver}_PerMaxDev',
                'std':f'BOSSFlats_{obsf.upper()}_{ver}_STD',
                'mad':f'BOSSFlats_{obsf.upper()}_{ver}_MAD',
                'count':f'BOSSFlats_{obsf.upper()}_{ver}_Count'}
    for fig in figs.keys():
        config = {'toImageButtonOptions': {'format': 'png','filename': fignames[fig],
                                           'scale': 6 }, # Multiply title/legend/axis/canvas sizes by this factor
                  'responsive': True}  # Ensure the figure is responsive

        fig_params = dict(full_html=False, default_height='900px', default_width='100%',
                          include_plotlyjs='cdn', config=config)
        figs[fig] = figs[fig].to_html(**fig_params)
    plotly_jinja_data = {"fig": figs['all'],
                         "fig_selct": figs['select'],
                         "fig_meta": meta.to_html(classes='scrollable-table',escape=False, index=False),
                         "fig_dev": figs['div'],
                         "fig_std": figs['std'],
                         "fig_mad": figs['mad'],
                         "fig_count": figs['count'],
                         "title":title, "name":name, "favicon":favicon}

    with open(filename, "w", encoding="utf-8") as output_file:
        with open(ptt.join(idlspec2d_dir,'templates','html','FieldLib_template.html')) as template_file:
            j2_template = Template(template_file.read())
            output_file.write(j2_template.render(plotly_jinja_data))

    

def csv_dump(Allflat, directory, version, fobs, meds_dic,  qbad_dic, beta_dic, confs_dic, confFiberids_dic):
    makedirs(ptt.join(directory,fobs,'csv'), exist_ok=True)
    fiberid = np.arange(1,501)
    exps=np.char.zfill(Allflat['EXP'].value.astype(str),8)
    rot=Allflat['ROT'].value
    alt=Allflat['ALT'].value
    obs=Allflat['OBS'].value
    mjd=Allflat['MJD'].value
    try: tai=Allflat['TAI'].value
    except: tai=np.full_like(Allflat['MJD'].value.tolist(), np.NaN)
    ccds = ['b1','r1'] if fobs == 'apo' else ['b2','r2']

    for i,row in enumerate(tqdm(meds_dic[ccds[0]], desc='CSV EXP', leave=False, position = 1)):
        flat_meta=pd.DataFrame()
        for ccd in ccds:
            flat_meta=pd.concat([flat_meta,pd.DataFrame({'FIBERID':fiberid, 'CONFFIBERID': confFiberids_dic[ccd][i],
                                                     'EXPID':[exps[i]]*500,
                                                     'SPECTROGRAPH':[ccd]*500,'CONFIGID':confs_dic[ccd][i],
                                                     'MED_FLAT_VALUE':meds_dic[ccd][i], 'ROT':[rot[i]]*500,
                                                     'ALT':[alt[i]]*500, 'OBS':np.array([obs[i]]*500),
                                                     'MJD':[mjd[i]]*500, 'TAI':[tai[i]]*500,
                                                     'BETA':beta_dic[ccd][i]})], ignore_index=False)
        flat_meta.to_csv(ptt.join(directory,fobs,'csv','summary_'+exps[i]+'.csv'),index=False)

def analysis(directory, version, mjd=None, noplot=False, obs='apo',
             lowFiber=.8, run='all', TraceIDs=False):
    Allflat=Table(fits.getdata(ptt.join(directory,'calibs','allflats.fits'),1))
    Allflat = Allflat[Allflat['OBS'] == obs]
    Allflat= tabunique(Allflat, keys='EXP')
    Allflat.sort('EXP')
    #Allflat = Allflat[:10]
    if run in ['all','med']:
        if mjd is None:
            outfile=ptt.join(directory,"FlatAnalysis-"+version+"_"+obs+".html")
        else:
            Allflat=Allflat[np.where(np.array(Allflat['MJD'].value.tolist()).astype(int) == int(mjd))]
            outfile=ptt.join(directory,"FlatAnalysis-"+version+"-mjd"+str(mjd)+"_"+obs+".html")
    if max(np.asarray(Allflat['MJD'].value).astype(int)) < 59550:
        plates = True
    else: plates = False
    meds_dic={}
    qbad_dic={}
    beta_dic={}
    assigns_dic={}
    confs_dic={}
    confFiberids_dic={}
    ccds = ['b1','r1'] if obs == 'apo' else ['b2','r2']
    if run in ['all','med','pos','csv','lowfiber']:
        for ccd in tqdm(ccds, desc='CCDs', position=0):
            meds=None
            qbad=[]
            confFiberids = None
            betas=None
            confs=None
            for i, exp in enumerate(tqdm(Allflat['EXP'].value.tolist(),desc='Exposure', position=1, leave=False)):
                
                if obs not in assigns_dic.keys():
                    assigns_dic[obs]=read_fiberAssignments(ptt.join(directory,'fiberAssignments',obs,'fiberAssignments.csv'))
                ff = ptt.join(directory,'calibs', obs, str(Allflat['MJD'].value.tolist()[i]),'spFlat-'+ccd+'-'+str(exp).zfill(8)+'.fits.gz')
                try:
                    data=fits.getdata(ff,0)
                    mid=data.shape[1]//2
                    data=data[:,mid-1000:mid+1000]
                    med=np.median(data,axis=1)
                except:
                    data=None
                    med = np.full(500,np.NaN)
                if data is not None:
                    hdr=fits.getheader(ff,0)
                    try:
                        confid = hdr['CONFID']
                    except:
                        confid = None
                    if (str(confid) == '-999') or (str(confid).strip().lower() == 'nan'):
                        confid = None
                    if confid is not None:
                        try:
                            confSummary = find_confSummary(confid, obs)
                        except:
                            print(ff)
                            raise
                    else:
                        confSummary = None
                else:
                    confid = None
                    confSummary  = None
                if confSummary is not None:
                    if len(confSummary.colnames) == 0:
                        confSummary = None
                if confSummary is not None:
                    beta = confSummary['beta'].value
                    confFiberid = confSummary['fiberId'].value
                else:
                    beta = np.full(500,np.NaN)
                    confFiberid = np.full(500,np.NaN)

                if betas is None:
                    meds = med
                    betas = beta
                    confFiberids = confFiberid
                    confs = np.array([confid]*len(beta))
                else:
                    meds=np.vstack([meds,med])
                    betas=np.vstack([betas,beta])
                    confFiberids=np.vstack([confFiberids,confFiberid])
                    confs=np.vstack([confs,np.array([confid]*len(beta))])

                if run in ['all','med','pos','lowfiber']:
                    ID_lowFiber(ff, ptt.join(directory, 'calibs', obs, 'lowFlagged', str(Allflat['MJD'].value.tolist()[i])) ,med, confFiberids, threshold=lowFiber)
                if run in ['lowfiber','csv']: continue
                if ptt.exists(ff) :
                    qbad.append(0)
                    if noplot is False:
                       mjdstr=str(Allflat['MJD'].value.tolist()[i])
                       makedirs(ptt.join(directory, 'calibs', obs, mjdstr),exist_ok=True)
                       makedirs(ptt.join(directory, 'calibs', obs, 'transmission', mjdstr),exist_ok=True)
                       if not ptt.exists(build_savename(ff,ptt.join(directory, 'calibs', obs, mjdstr))):
                            plot_flat(ff,ptt.join(directory, 'calibs', obs, mjdstr))
                            raw_flat=get_raw(ff,obs, mjdstr)
                            plot_raw(raw_flat,ptt.join(directory, 'calibs', obs, mjdstr))
                            plot_thruput_v_sextant(ff,str(Allflat['MJD'].value.tolist()[i]), ptt.join(directory, 'calibs', obs, 'transmission', mjdstr), assigns_dic[obs])
                            plot_thruput_v_sextant(ff,str(Allflat['MJD'].value.tolist()[i]), ptt.join(directory, 'calibs', obs, mjdstr), assigns_dic[obs])
                else: qbad.append(1)
            if run not in ['csv']:
                ID_lowFiber('merged_'+ccd, ptt.join(directory, 'calibs',obs, 'lowFlagged') ,meds, confFiberids, threshold=lowFiber, clobber=True)
            meds_dic[ccd]=meds
            qbad_dic[ccd]=np.asarray(qbad)
            beta_dic[ccd]=betas
            confs_dic[ccd]=confs
            confFiberids_dic[ccd]=confFiberids
            if run  in ['lowfiber','csv']: continue
        if run == 'lowfiber': return

        if run in ['all', 'med']:
            tqdm.write('Producing med_plot')
            
            tabs=med_plot(Allflat,ptt.join(directory,'calibs'),
                          ptt.join(directory,f"FlatAnalysis-{version}_{obs}.html"),
                          version, obs, meds_dic[ccds[0]],
                          meds_dic[ccds[1]], qbad_dic[ccds[0]],
                          qbad_dic[ccds[1]], confFiberids_dic[ccds[0]])
        if run in ['all', 'pos']:
            position_plots(Allflat, ptt.join(directory,'calibs'), version, obs,
                            meds_dic[ccds[0]], meds_dic[ccds[1]],
                            qbad_dic[ccds[0]], qbad_dic[ccds[1]],
                            beta_dic[ccds[0]], beta_dic[ccds[1]],
                            confFiberids_dic[ccds[0]], plates=plates)
        if run in ['all','csv']:
            csv_dump(Allflat, ptt.join(directory,'calibs'), version, obs, meds_dic,  qbad_dic, beta_dic, confs_dic,confFiberids_dic)
    if run in ['timeSeries', 'all']:
        plot_thruput_timeseries(ptt.join(directory,'calibs', obs,'csv'), Traceid=TraceIDs)
