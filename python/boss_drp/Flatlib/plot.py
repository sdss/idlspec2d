#!/usr/bin/env python3
import os.path as ptt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from os import getenv, makedirs, sep, rename
from astropy.io import fits
import argparse
import glob
import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

def plot_flat(filename, save_dir):
    loc=get_mdj_obs(filename)
    outname, shortname = build_savename(filename,save_dir,return_short=True)
    try:
        flat_arr=fits.getdata(filename,0)
    except:
        print('error with '+filename)
        flat_arr = [[np.NaN,np.NaN],[np.NaN,np.NaN]]
    plt.figure(figsize=(12, 2), dpi=100)
    plt.imshow(flat_arr,cmap='gray')
    plt.xlabel('pixel')
    plt.ylabel('Fiber')
    plt.title('MJD:'+str(loc['mjd'])+' '+shortname)
    plt.savefig(outname,bbox_inches='tight')
    plt.close('all')
    return
        
def build_savename(filename, save_dir,return_short=False, mod=''):
    makedirs(save_dir, exist_ok=True)
    filename=ptt.basename(filename)
    filename=ptt.splitext(ptt.splitext(filename)[0])[0]
    if return_short is True:
        return(ptt.join(save_dir, filename+mod+'.png'),filename)
    return(ptt.join(save_dir, filename+mod+'.png'))


def read_fiberAssignments(filename):
    if not ptt.exists(filename):
        makedirs(ptt.dirname(filename), exist_ok=True)
        if 'lco' in filename:
            url = 'https://raw.githubusercontent.com/sdss/fps_calibrations/main/lco/wok_calibs/duPontFlatCMM/fiberAssignments.csv'
        else:
            url = 'https://raw.githubusercontent.com/sdss/fps_calibrations/main/apo/wok_calibs/sloanFlatCMM/fiberAssignments.csv'
        response = requests.get(url)
        if response.status_code == 200:
            # Write the content to a local file
            with open(filename, 'wb') as f:
                f.write(response.content)
        else:
            raise Exception(f'Failed to download file: {response.status_code}')
    
    fiberAssign=pd.read_csv(filename,index_col=0)
    fiberAssign=fiberAssign[np.isfinite(fiberAssign.BOSSFiber)]
    fiberAssign=fiberAssign.sort_values(by=['BOSSFiber'])
    return(fiberAssign)

def plot_thruput_v_sextant(filename,mjd, save_dir, fiberAssignments):
    outname, shortname = build_savename(filename, save_dir, return_short=True, mod='_trans')
    def reject_outliers(data, m=2):
        return data[abs(data - np.mean(data)) < m * np.std(data)]

    try:
        data=fits.getdata(filename,0)
    except:
        print('error with '+filename)
        data = [[np.NaN,np.NaN],[np.NaN,np.NaN]]
    mean=np.mean(data,axis=1)
    arr0=mean
    arr=(reject_outliers(mean,m=3))
    filtdata=data
    while np.abs(np.mean(arr)-np.mean(arr0)) > .0001:
        arr0=arr
        arr=reject_outliers(arr,m=3)
    trans=mean/np.mean(arr)

    ax=plt.subplot()
    for i in range(1, 7):
        filt=np.where(fiberAssignments.Sextant.values == i)[0]
        im=ax.scatter(fiberAssignments.Sextant.values[filt], trans[filt], marker='.',
                        lw=.1,ec='k',c=np.arange(len(filt)), cmap='CMRmap')
    ax.set_xlabel('Sextant')
    ax.set_ylabel('Relative Transmission (to mean)')
    plt.title('MJD:'+str(mjd)+' '+shortname)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar=plt.colorbar(im, cax=cax)
    cbar.set_label('Fiber within Sextant')
    plt.savefig(outname, bbox_inches='tight')
    plt.close('all')
    return

def uni(lst):
    arr = np.array(lst, dtype=object)  # Keep NaNs as objects
    unique_vals = list(set(arr[~np.isnan(arr.astype(float))]))  # Unique non-NaNs
    if np.any(np.isnan(arr.astype(float))):  # Check if NaN exists
        unique_vals.append(np.nan)
    return unique_vals

def plot_thruput_timeseries(csv_dir, Traceid=False):
    files = glob.glob(ptt.join(csv_dir,'summary_*.csv'))
    ts_r = None
    ts_b = None
    ids_b = None
    ids_r = None
    for file in tqdm(files, desc='Reading CSVs', leave=False, position=1):
        df = pd.read_csv(file)
        df_b = df[df['SPECTROGRAPH'].isin(['b1','b2'])]
        df_r = df[df['SPECTROGRAPH'].isin(['r1','r2'])]

        if ts_r is None:
            ids_b = df_b['CONFFIBERID'].tolist() if not Traceid else df_b['FIBERID'].tolist()
            ids_r = df_r['CONFFIBERID'].tolist() if not Traceid else df_r['FIBERID'].tolist()

            ts_b = pd.DataFrame([df_b['MED_FLAT_VALUE'].tolist()], index=[df.iloc[0]['MJD']], columns = ids_b)
            ts_r = pd.DataFrame([df_r['MED_FLAT_VALUE'].tolist()], index=[df.iloc[0]['MJD']], columns = ids_r)

        else:
            ids_b = df_b['CONFFIBERID'].tolist() if not Traceid else df_b['FIBERID'].tolist()
            ids_r = df_r['CONFFIBERID'].tolist() if not Traceid else df_r['FIBERID'].tolist()
            mjd_u = df.iloc[0]['MJD']
            if (len(uni(ids_b)) != len(ids_b)) or (len(uni(ids_r)) != len(ids_r)): continue
            if mjd_u not in ts_b.index:
                try:
                    ts_b.loc[mjd_u] = dict(zip(ids_b, df_b['MED_FLAT_VALUE'].tolist()))
                except Exception as e:
                    tqdm.write(f'Error with {file}: {e}')
            else:
                temp = ts_b.loc[mjd_u].copy()
                new_data = dict(zip(ids_b, df_b['MED_FLAT_VALUE'].tolist()))
                for key, value in new_data.items():
                    if key not in ts_b.columns:
                        if np.isnan(key):
                            continue

                        tqdm.write(f'adding {key} - {file}')
                        ts_b = ts_b.copy()
                        ts_b[key] = np.nan  # Ensure the column exists
                    temp[key] = value
                ts_b.loc[mjd_u] = temp
            if mjd_u not in ts_r.index:
                try:
                    ts_r.loc[mjd_u] = dict(zip(ids_r, df_r['MED_FLAT_VALUE'].tolist()))
                except Exception as e:
                    tqdm.write(f'Error with {file}: {e}')
            else:
                temp = ts_r.loc[mjd_u].copy()
                new_data = dict(zip(ids_r, df_r['MED_FLAT_VALUE'].tolist()))
                for key, value in new_data.items():
                    if key not in ts_r.columns:
                        if np.isnan(key):
                            continue
                        tqdm.write(f'adding {key} - {file}')
                        ts_r = ts_r.copy()
                        ts_r[key] = np.nan  # Ensure the column exists
                    temp[key] = value
                ts_r.loc[mjd_u] = temp

    ts_b = ts_b.sort_index()
    ts_r = ts_r.sort_index()
    makedirs(ptt.join(ptt.dirname(csv_dir),'timeseries'), exist_ok=True)
    if Traceid:
        flag = 'Trace ' #fiberid on chip
    else:
        flag = 'Fiber ' #Fiberid on slit
    cols = ts_b.columns.tolist()+ts_r.columns.tolist()
    cols = np.sort(np.unique(np.asarray(cols))).astype(int).tolist()
    for i in range(len(cols)):
        if((i % 10) == 0):
            fig, axs= plt.subplots(2,1,figsize=(10,5), dpi =100)
            s =i
        axs[0].plot(ts_b.index.tolist(), ts_b[cols[i]], marker='.', label=f'{flag}{cols[i]}')
        axs[1].plot(ts_r.index.tolist(), ts_r[cols[i]], marker='.', label=f'{flag}{cols[i]}')
        axs[0].set_ylabel('Blue Relative \nTransmission (to mean)')
        axs[1].set_ylabel('Red Relative \nTransmission (to mean)')
        axs[1].set_xlabel('MJD')
        axs[0].set_ylim([-.05,1.2])
        axs[1].set_ylim([-.05,1.2])
        axs[0].legend(ncols=10, fontsize=7, framealpha=.2,bbox_to_anchor=(0.5, 1.15),loc='upper center')#,bbox_transform=axs[0].transAxes,)
        axs[1].legend(ncols=10, fontsize=7, framealpha=.2,bbox_to_anchor=(0.5, 1.15),loc='upper center')#,bbox_transform=axs[0].transAxes,)
#        axs[0].legend(ncols=10, fontsize=8, framealpha=.2,bbox_to_anchor=(0.5, 1.05),loc='center',bbox_transform=axs[0].transAxes,)
#        axs[1].legend(ncols=10, fontsize=8, framealpha=.2,bbox_to_anchor=(0.5, 1.05),loc='center',bbox_transform=axs[0].transAxes,)
        if((i % 10) == 9) or (i==len(cols)-1):
            fig.tight_layout()
            plt.subplots_adjust()#top=.93)  # Adjust top to fit the legends
            if Traceid:
                plt.savefig(ptt.join(ptt.dirname(csv_dir),'timeseries',
                            f'timeseries_Tracefiber{str(s+1).zfill(3)}-{str(i+1).zfill(3)}.png'),bbox_inches='tight')
            else:
                plt.savefig(ptt.join(ptt.dirname(csv_dir),'timeseries',
                            f'timeseries_fiber{str(s+1).zfill(3)}-{str(i+1).zfill(3)}.png'),bbox_inches='tight')
            plt.close('all')


    if Traceid:
        def key(x):
            x = ptt.basename(x)
            x = x.replace('timeseries_Tracefiber','').split('-')[0]
            return(int(x))
        build_preview(ptt.dirname(csv_dir),'timeseries_Tracefiber', deep=1, nper=2,pattern='timeseries_Tracefiber*png', key=key)
        with open(ptt.join(ptt.dirname(csv_dir),'timeseries_Tracefiber_data_blue.html'),'w') as f: f.write(ts_b.to_html())
        with open(ptt.join(ptt.dirname(csv_dir),'timeseries_Tracefiber_data_red.html'),'w') as f: f.write(ts_r.to_html())
        ts_r.to_csv(ptt.join(ptt.dirname(csv_dir),'timeseries_Tracefiber_data_red.csv'), index=False)
        ts_b.to_csv(ptt.join(ptt.dirname(csv_dir),'timeseries_Tracefiber_data_blue.csv'), index=False)
    else:
        def key(x):
            x = ptt.basename(x)
            x = x.replace('timeseries_fiber','').split('-')[0]
            return(int(x))
        build_preview(ptt.dirname(csv_dir),'timeseries', deep=1, nper=2,pattern='timeseries_fiber*png', key=key)
        with open(ptt.join(ptt.dirname(csv_dir),'timeseries_data_red.html'),'w') as f: f.write(ts_r.to_html())
        with open(ptt.join(ptt.dirname(csv_dir),'timeseries_data_blue.html'),'w') as f: f.write(ts_b.to_html())
        ts_r.to_csv(ptt.join(ptt.dirname(csv_dir),'timeseries_data_red.csv'), index=False)
        ts_b.to_csv(ptt.join(ptt.dirname(csv_dir),'timeseries_data_blue.csv'), index=False)


def plot_raw(filename, save_dir):
    loc=get_mdj_obs(filename)
    img_arr=fits.getdata(filename,0)
    filename=ptt.basename(filename)
    filename=ptt.splitext(ptt.splitext(filename)[0])[0]
    plt.figure(figsize=(6, 6), dpi=100)
    plt.imshow(img_arr,cmap='gray')
    plt.xlabel('pixel')
    plt.ylabel('pixel')
    plt.title('MJD:'+str(loc['mjd'])+' '+filename)
    plt.savefig(ptt.join(save_dir, filename+'.png'),bbox_inches='tight')
    plt.close('all')
    return
    
def get_mdj_obs(filename):
    mjd=ptt.basename(ptt.dirname(ptt.normpath(filename)))
    obs=ptt.basename(ptt.dirname(ptt.dirname(ptt.normpath(filename))))
    return({'mjd':mjd,'obs':obs})

def get_raw(filename, obs, mjd):
    filename=ptt.basename(filename)
    filename=filename.replace('spFlat-','sdR-')
    filename=filename.replace('.fits.gz','.fit.gz')
    if obs.lower() == 'apo': env='BOSS_SPECTRO_DATA_N'
    else:  env='BOSS_SPECTRO_DATA_S'
    filename=ptt.join(getenv(env),str(mjd),filename)
    return(filename)

def build_preview(directory, name, deep=1, nper=2,pattern='spFlat-??-????????.png', key = None):
    sub = sep.join(['*']*deep)
    figs = glob.glob(ptt.join(directory,sub,pattern))
    if key is not None:
        figs = sorted(figs, key=key)

    with open(ptt.join(directory, 'tmp-index.html'), 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>'+'\n')
        f.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'+'\n')
        f.write('<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">'+'\n')
        f.write('<head>'+'\n')
        f.write('<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />'+'\n')
        f.write('<title>'+name+'</title>'+'\n')
        f.write('<style type="text/css">'+'\n')
        f.write('body { background: #111; }'+'\n')
        f.write('td { color: gray; }'+'\n')
        f.write('</style>'+'\n')
        f.write('</head>'+'\n')
        f.write('<body>'+'\n')
        f.write('<table border="0" cellspacing="5">'+'\n')

 

        for i,ff in enumerate(figs):
            ff = ptt.relpath(ff,start = directory)
            f.write(f'<td>{ff}<br><a href="{ff}"><img src="{ff}" ></a></td>\n')
            if(((i % nper) == nper-1) or (i == len(figs)-1)): f.write('</tr>'+'\n')
        f.write('</table>'+'\n')
        f.write('</body>'+'\n')
        f.write('</html>'+'\n')
    rename(ptt.join(directory, 'tmp-index.html'), ptt.join(directory, f'{name}.html'))

def plot(dir_, flat, savedir ='.', assigns=None):
    if assigns is None:
        assigns = {}
        for obs in ['apo','lco']:
            assigns[obs]=read_fiberAssignments(ptt.join(dir_,'fiberAssignments',obs,'fiberAssignments.csv'))
    loc=get_mdj_obs(flat)
    makedirs(ptt.join(savedir,loc['mjd']), exist_ok=True)
    makedirs(ptt.join(savedir,'transmission',loc['mjd']), exist_ok=True)
    
    plot_flat(flat,ptt.join(savedir,loc['mjd']))
    build_preview(savedir, 'flats', deep=1)
    
    raw_flat = get_raw(flat, loc['obs'], loc['mjd'])
    plot_raw(raw_flat, ptt.join(savedir,loc['mjd']))
    build_preview(savedir, 'raw', deep=1, pattern='sdR*.png', nper=4)
    
    plot_thruput_v_sextant(flat,loc['mjd'], ptt.join(savedir, 'transmission',loc['mjd']), assigns[loc['obs']])
    build_preview(ptt.join(savedir, 'transmission'), 'transmission', deep=1,pattern='spFlat*_trans.png',nper=2)
    return
