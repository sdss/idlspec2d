#!/usr/bin/env python3

import argparse
import sys
import os.path as ptt
import logging
import numpy as np
from pydl.pydlutils.yanny import yanny, read_table_yanny
from astropy.io import fits
from astropy.table import Table, Column, unique
from os import getenv, remove, makedirs
from glob import glob
import time
import datetime
import astropy.time
import os
import re
import time
from pydl.pydlutils import sdss
import mmap
from merge_dm import merge_dm
from splog import Splog
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
splog = Splog()


# this version does not 
# - remove partial epochs
# - allow purge of partial
# - check for aborted combines

oplimits = None
oplimit_filename=''
chunkdata = None
pulic_plate_data = None
spPlatelistMessage =False
try:
    sdss.set_maskbits(maskbits_file=getenv("IDLUTILS_DIR")+"/data/sdss/sdssMaskbits.par")
except:
    splog.log('Environmental Varable IDLUTILS_DIR must be set')
    exit()


def plot_sky(topdir, flist, flist_file):
    idx = np.where(np.char.strip(flist['STATUS1D'].data) == 'Done')[0]

    RA   = flist['RACEN'].data[idx]
    DEC  = flist['DECCEN'].data[idx]
    prog = np.char.upper(flist[idx]['PROGRAMNAME'].data)
    fcad = np.char.lower(flist[idx]['FIELD_CADENCE'].data)
    fsur = np.char.lower(flist[idx]['SURVEY'].data)
    status = np.char.lower(flist[idx]['STATUS1D'].data)

# plot the RA/DEC in an area-preserving projection
# convert coordinates to degrees
    RA *= np.pi / 180
    DEC *= np.pi / 180
    phi = np.linspace(0, 2.*np.pi, 36)  #36 points
    r = np.radians(1.5)
    C0=C1=C2=C3=C4=C5=C6=C7=0

    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.0))
    ax = plt.axes(projection='mollweide')
    plt.grid(True)
    plt.title('SDSS plate/field locations')
    for i in range(0, len(RA)):
        
        if RA[i] < np.pi:
            x = RA[i] + r*np.cos(phi)
        else:
            x = RA[i] + r*np.cos(phi)-2*np.pi
        y = DEC[i] + r*np.sin(phi)
        if 'dark' in fcad[i].lower():
            label = 'FPS DARK'  if C0 == 0 else None
            ax.plot(x, y, color = 'blue', label = label, alpha=.5)
            C0 = 1
        elif 'bright' in fcad[i]:
            label = 'FPS BRIGHT'  if C1 == 0 else None
            ax.plot(x, y, color = 'orange', label = label, alpha=.5)
            C1 = 1
        elif 'plate' in fcad[i]:
            label = 'PLATE'  if C2 == 0 else None
            ax.plot(x, y, color = 'green', label = label, alpha=.5)
            C2 = 1
        else:
#            pass
            label = 'MANUAL'  if C3 == 0 else None
            ax.plot(x, y, color = 'saddlebrown', label = label, alpha=.5)
            C3 = 1

    plt.legend(loc=1,fontsize=10)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab)
    fig.tight_layout()
    plt.savefig(ptt.join(topdir,'SDSSVc_s.png'),dpi=50,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSVc.png'),dpi=500,bbox_inches='tight')
    plt.close()

    C0=C1=C2=C3=C4=C5=C6=C7=C8=C9=0

####################################################################################
    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.0))
    ax = plt.axes(projection='mollweide')
    plt.grid(True)
    plt.title('SDSS-V plate/field locations')
    for i in range(0, len(RA)):
        if RA[i] < np.pi:
            x = RA[i] + r*np.cos(phi)
        else:
            x = RA[i] + r*np.cos(phi)-2*np.pi
        y = DEC[i] + r*np.sin(phi)
        if 'RM' in prog[i]:
            label = 'RM'  if C0 == 0 else None
            ax.plot(x, y, color = 'blue', label = label)
            C0 = 1
        elif 'MWM' in prog[i]:
            label = 'MWM'  if C1 == 0 else None
            ax.plot(x, y, color = 'red', label = label)
            C1 = 1
        elif 'AQMES-Medium'.upper() in prog[i]:
            label = 'AQMES-Medium'  if C2 == 0 else None
            ax.plot(x, y, color = 'pink', label = label)
            C2 = 1
        elif 'AQMES-Wide'.upper() in prog[i]:
            label = 'AQMES-Wide'  if C3 == 0 else None
            ax.plot(x, y, color = 'orange', label = label)
            C3 = 1
        elif 'AQMES-Bonus'.upper() in prog[i]:
            label = 'AQMES-Bonus'  if C4 == 0 else None
            ax.plot(x, y, color = 'magenta', label = label)
            C4 = 1
        elif 'eFEDS'.upper() in prog[i]:
            label = 'eFEDS' if C5 == 0 else None
            ax.plot(x, y, color = 'green', label = label)
            C5 = 1
        elif 'OFFSET'.upper() in prog[i]:
            label = 'OFFSET'  if C6 == 0 else None
            ax.plot(x, y, color = 'black', label = label)
            C6 = 1
        elif 'mwm-bhm' in fsur[i]:
            label = 'MWM-BHM FPS'  if C7 == 0 else None
            ax.plot(x, y, color = 'navy', label = label)
            C7 = 1
        elif 'bhm-mwm' in fsur[i]:
            label = 'BHM-MWM FPS'  if C8 == 0 else None
            ax.plot(x, y, color = 'lightcoral', label = label)
            C8 = 1
        else:
            ax.plot(x, y, color = 'saddlebrown')
            C9 = 1
    plt.legend(loc=1,fontsize=10)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab)
    fig.tight_layout()
    plt.savefig(ptt.join(topdir,'SDSSV_s.png'),dpi=50,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSV.png'),dpi=500,bbox_inches='tight')
    plt.close()

####################################################################################
    file_c = flist_file.replace('fieldlist','spAll')+'.gz'
    if os.path.exists(file_c):
        hdu_list = fits.open(file_c)
        table_hdu = hdu_list[1]
        table_data = table_hdu.data
    else: table_data=Table(names=('RACAT','DECCAT', 'PROGRAMNAME','NSPECOBS'))
    RA1=table_data.field('RACAT')
    DEC1=table_data.field('DECCAT')
    Nobs=table_data.field('NSPECOBS')

    if len(RA1) == 0:
        RA1 =np.asarray([np.NaN])
        DEC1=np.asarray([np.NaN])
        Nobs=np.zeros_like(DEC1)
    RA1 *= np.pi / 180
    DEC1 *= np.pi / 180
    for i in range(0, len(RA1)):
        if RA1[i] >= np.pi:
            RA1[i]=RA1[i]-2*np.pi

    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.0))
    ax = plt.axes(projection='mollweide')
    plt.scatter(RA1, DEC1, s=2, c=Nobs, cmap=plt.cm.jet, edgecolors='none', linewidths=0)
    plt.grid(True)
    plt.title('SDSS Observed Targets')
    cb = plt.colorbar(cax=plt.axes([0.05, 0.1, 0.9, 0.05]),
                      orientation='horizontal')
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab)
    cb.set_label('Nobs')
    plt.savefig(ptt.join(topdir,'SDSSV2_s.png'),dpi=50,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSV2.png'),dpi=500,bbox_inches='tight')
    plt.close()

####################################################################################
    allpointings = Table()
    idx = np.where(np.char.strip(flist['STATUS1D'].data) == 'Done')[0]
    if len(allpointings) == 0:
        allpointings['RA'] = [np.NaN]
        allpointings['DEC'] = [np.NaN]
        allpointings['Nexp'] = [np.NaN]
    else:
        allpointings['RA']  = flist['RACEN'].data[idx]
        allpointings['DEC'] = flist['DECCEN'].data[idx]
        allpointings['Nexp'] = flist['NEXP'].data[idx]
    
        
    pointings = unique(allpointings, keys=['RA','DEC'])
    pointings.add_column(0, name='Nobs')
    pointings['Nexp'] = 0

    pointings.add_column(0.0, name='x')
    pointings.add_column(0.0, name='y')
    for row in pointings:
        idx = np.where((allpointings['RA'].data == row['RA']) &
                       (allpointings['DEC'].data == row['DEC']))[0]
        row['Nobs'] = len(idx)
        row['Nexp'] = np.sum(allpointings[idx]['Nexp'].data)
        if row['RA']*np.pi / 180 < np.pi:
            row['x'] = row['RA']*np.pi / 180
        else:
            row['x'] = row['RA']*np.pi / 180 - 2*np.pi
        row['y'] = row['DEC']*np.pi / 180
                
    
    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.0))
    ax = plt.axes(projection='mollweide')
    plt.scatter(pointings['x'].data, pointings['y'].data, s=30,
                c=pointings['Nobs'].data, marker= 'x',
                cmap=plt.cm.jet,alpha=.5)

    plt.grid(True)
    plt.title('SDSS field locations')
    cb = plt.colorbar(cax=plt.axes([0.05, 0.1, 0.9, 0.05]),
                      orientation='horizontal')
        
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab)
    cb.set_label('Number of Field-MJD ')
    plt.savefig(ptt.join(topdir,'SDSSV3.png'),dpi=500,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSV3_s.png'),dpi=50,bbox_inches='tight')
    plt.close()


def field_to_string(field):
    return(str(int(field)).zfill(6))

def getquality(row,basehtml,epoch=False, dereddened_sn2=False, rawsn2=False):
    sfield = field_to_string(row['FIELD'])
    if epoch is True:
        row['PLOTSN'] = '<a href="'+basehtml+'/'+row['RUN2D']+'/'+sfield+'/epoch/spSN2d-'+sfield+'-'+str(row['MJD'])+'.ps">SNPLOT</a>'
    else:
        row['PLOTSN'] = '<a href="'+basehtml+'/'+row['RUN2D']+'/'+sfield+'/spSN2d-'+sfield+'-'+str(row['MJD'])+'.ps">SNPLOT</a>'
    row['FIELDQUALITY'] = ''
    
    
    NEXP_B1 = row['NEXP_B1'] if 'NEXP_B1' in row else 0
    NEXP_R1 = row['NEXP_R1'] if 'NEXP_R1' in row else 0
    NEXP_B2 = row['NEXP_B2'] if 'NEXP_B2' in row else 0
    NEXP_R2 = row['NEXP_R2'] if 'NEXP_R2' in row else 0
    
    SN2_G1 = row['SN2_G1'] if 'SN2_G1' in row else 0
    SN2_I1 = row['SN2_I1'] if 'SN2_I1' in row else 0
    SN2_G2 = row['SN2_G2'] if 'SN2_G2' in row else 0
    SN2_I2 = row['SN2_I2'] if 'SN2_I2' in row else 0
 
    DERED_SN2_G1 = row['DERED_SN2_G1'] if 'DERED_SN2_G1' in row else 0
    DERED_SN2_I1 = row['DERED_SN2_I1'] if 'DERED_SN2_I1' in row else 0
    DERED_SN2_G2 = row['DERED_SN2_G2'] if 'DERED_SN2_G2' in row else 0
    DERED_SN2_I2 = row['DERED_SN2_I2'] if 'DERED_SN2_I2' in row else 0
 
    nexps = np.array([NEXP_B1,NEXP_R1,NEXP_B2,NEXP_R2])
    valid = np.where(nexps !=0)[0]
    valid_spb = np.where([NEXP_B1,NEXP_B2])[0]
    valid_spr = np.where([NEXP_R1,NEXP_R2])[0]
    nexp_max  = 0
    nexp_min  = 0
    
    if len(valid) != 0:
        nexp_min = min(nexps[valid])
        nexp_max = max(nexps[valid])
        fieldsn2 = np.array([SN2_G1,SN2_I1,SN2_G2,SN2_I2])
        row['FIELDSN2'] = min(fieldsn2[valid])
        
        if dereddened_sn2:
            deredsn2 = np.array([DERED_SN2_G1,DERED_SN2_I1,DERED_SN2_G2,DERED_SN2_I2])
            row['DEREDSN2'] = min(deredsn2[valid])

        
        sn2b = np.array([SN2_G1,SN2_G2])
        sn2r = np.array([SN2_I1,SN2_I2])
        if dereddened_sn2:
            sn2b = np.array([DERED_SN2_G1,DERED_SN2_G2])
            sn2r = np.array([DERED_SN2_I1,DERED_SN2_I2])
        
        min_sn2_b = min(sn2b[valid_spb]) if len(valid_spb) != 0 else 0
        min_sn2_r = min(sn2r[valid_spr]) if len(valid_spr) != 0 else 0
    iqual = 2
    
    prog = row['PROGRAMNAME'].strip() if 'PROGRAMNAME' in row else ''
    mjd  = row['MJD']

    is_elg_plate = True if prog.upper()  in ['ELG_NGC','ELG_SGC'] else False

    if int(row['FIELD']) < 15000:
        #--- JEB 2018-05-23: if elg plate, plate is 'good' no matter what SN2
        if not is_elg_plate:
            #--- JEB 2018-05-23: new thresholds after 2017-10-03
            if (int(mjd) > 58029):
                if ((min_sn2_b < 8.0) or (min_sn2_r < 18.0)):
                    iqual = min(iqual,0)
                if ((min_sn2_b < 10.0) or (min_sn2_r < 22.0)):
                    iqual = min(iqual,0)
    if row['FBADPIX'] > 0.10:
        iqual = iqual < 0
    # For reductions before v5_1, NEXP_MIN and NEXP_MAX are always zero
    if (nexp_max > 0):
        if int(row['FIELD']) < 16000:
            if (nexp_min < 3):
                iqual = min(iqual,0)
    min_sn2_b_scaled=min_sn2_b/(row['EXPTIME']/3600)
    min_sn2_r_scaled=min_sn2_r/(row['EXPTIME']/3600)
    if (((min_sn2_b < 10.) or (min_sn2_r < 22.)) &
        ((min_sn2_b_scaled < 10.0) or (min_sn2_r_scaled < 22.0))):
            iqual = min(iqual,0)
    qualstring = ['bad', 'marginal', 'good']
    if row['FIELDQUALITY'] == '':
        row['FIELDQUALITY'] = qualstring[iqual]

    return(row)

def getoutputs(row,basehtml,epoch=False):
    sfield = field_to_string(row['FIELD'])
    if epoch is True:
        row['PLOTS'] = '<a href="'+basehtml+'/images/'+row['RUN2D']+'/epoch/'+row['RUN1D']+'/'+sfield+'-'+str(row['MJD'])+'/">PLOTS</a>'
        row['DATA']  = '<a href="'+basehtml+'/'+row['RUN2D']+'/epoch/spectra/full/'+sfield+'/'+str(row['MJD'])+'/">DATA</a>'
    else:
        row['PLOTS'] = '<a href="'+basehtml+'/images/'+row['RUN2D']+'/'+row['RUN1D']+'/'+sfield+'-'+str(row['MJD'])+'/">PLOTS</a>'
        row['DATA']  = '<a href="'+basehtml+'/'+row['RUN2D']+'/spectra/full/'+sfield+'/'+str(row['MJD'])+'/">DATA</a>'

    return(row)
    
def get_chunkinfo(row):
    global chunkdata
    if chunkdata is None:
        chunkfile = ptt.join(getenv('PLATELIST_DIR'), 'platePlans.par')
        try:
            chunkdata = Table(yanny(chunkfile)['PLATEPLANS'])
        except:
            splog.log('Empty or missing platePlans.par file')
            chunkdata = None
            return(row)
    cinfo = chunkdata[np.where(chunkdata['plateid'] == int(row['FIELD']))[0]][0]
    row['SURVEY'] = cinfo['survey']
    row['PROGRAMNAME'] = cinfo['programname']
    row['CHUNK'] = cinfo['chunk']
    row['FILEID'] = cinfo['tileid']
    row['DESIGNS'] = str(cinfo['designid'])
    row['RACEN'] = cinfo['raCen']
    row['DECCEN'] = cinfo['decCen']
    row['EPOCH'] = cinfo['epoch']
    row['FIELD_CADENCE'] = 'Plates'
    row['CHUNKHTML'] = '<a href="https://platedesign.sdss.org/runs/'+cinfo['chunk']+'/'+cinfo['chunk']+'.html">'+cinfo['chunk']+'</a>'
    return(row)

def wwhere(array, value):
    value = value.replace('*','[\w]*')
    r = re.compile(value, re.IGNORECASE)
    ret = np.full(len(array), False)
    idx = [i for i, x in enumerate(array) if r.search(x)]
    ret[idx] = True
    return(ret)

def read_spec1d(row, path, fieldfile, epoch=False):
    spZfile     = ptt.join(path, row['RUN1D'], 'spZbest-'+row['FIELD']+'-'+row['MJD']+'.fits')
    spDiag1dlog = ptt.join(path, row['RUN1D'], 'spDiag1d-'+row['FIELD']+'-'+row['MJD']+'.log')
    if ptt.exists(spZfile):
        strplt = field_to_string(row['FIELD'])
        strmjd = str(row['MJD']).strip()
        run2d = row['RUN2D']
        run1d = row['RUN1D']
        
        try:
            zans = Table(fits.getdata(spZfile, 1))
            plug = Table(fits.getdata(fieldfile,5))
        except:
            zans = None
            plug = None
            pass
        if zans is not None:

            objclass = np.char.strip(np.char.upper(zans['CLASS'].data))
            objtyp = np.char.strip(np.char.upper(plug['OBJTYPE'].data))
            # Use the ZWARNING flag if it exists to identify SKY or UNKNOWN.
            if 'ZWARNING' in zans.columns:
                zwarning = zans['ZWARNING'].data
            else:
                zwarning = np.zeros(len(zans), dtype=int)
            qsky = (zwarning & 1) != 0
            row['N_GALAXY']  = len(np.where((objclass == 'GALAXY')   & (zwarning == 0) & (objtyp != 'SPECTROPHOTO_STD') & (objtyp != 'SKY'))[0])
            row['N_QSO']     = len(np.where((objclass == 'QSO')      & (zwarning == 0) & (objtyp != 'SPECTROPHOTO_STD') & (objtyp != 'SKY'))[0])
            row['N_STAR']    = len(np.where((objclass == 'STAR')     & (zwarning == 0) & (objtyp != 'SPECTROPHOTO_STD') & (objtyp != 'SKY'))[0])
            row['N_STD']     = len(np.where((objtyp   == 'SPECTROPHOTO_STD'))[0])
            row['N_UNKNOWN'] = len(np.where(((objclass == 'UNKNOWN') & (objtyp != 'SPECTROPHOTO_STD') & (objtyp != 'SKY')) | ((zwarning !=0) & (qsky == 0)))[0])
            row['N_SKY']     = len(np.where((objclass == 'SKY') | (qsky == 1))[0])
            row['STATUS1D']  = 'Done'


                
            
            nobj = len(zans)
            target = np.full(nobj,'', dtype=object)
            if int(row['FIELD']) < 15000:
                if 'PRIMTARGET' in np.char.upper(plug.columns):
                    for i in range(nobj):
                        target[i]+=sdss.sdss_flagname('TARGET', plug['PRIMTARGET'], concat = True)+' '
                if 'BOSS_TARGET1' in np.char.upper(plug.columns):
                    for i in range(nobj):
                        target[i]+=sdss.sdss_flagname('BOSS_TARGET1', plug['BOSS_TARGET1'], concat = True)+' '
                if 'EBOSS_TARGET0' in np.char.upper(plug.columns):
                    for i in range(nobj):
                        target[i]+=sdss.sdss_flagname('EBOSS_TARGET0', plug['EBOSS_TARGET0'], concat = True)+' '
                if 'EBOSS_TARGET1' in np.char.upper(plug.columns):
                    for i in range(nobj):
                        target[i]+=sdss.sdss_flagname('EBOSS_TARGET1', plug['EBOSS_TARGET1'], concat = True)+' '
                if 'EBOSS_TARGET2' in np.char.upper(plug.columns):
                    for i in range(nobj):
                        target[i]+=sdss.sdss_flagname('EBOSS_TARGET2', plug['EBOSS_TARGET2'], concat = True)+' '
#            elif int(row['FIELD']) < 16000:
#                if 'SDSSV_BOSS_TARGET0' in np.char.upper(plug.columns):
#                    for i in range(nobj):
#                        target[i]+=sdss.sdss_flagname('SDSSV_BOSS_TARGET0', plug['SDSSV_BOSS_TARGET0'], concat = True)+' '
            for i in range(nobj):
                if target[i] == '': target[i] = zans['CLASS'][i]
            
            # Objects which shouldn't count against the success statistics
            bad_fiber = sdss.sdss_flagval('ZWARNING', 'LITTLE_COVERAGE') | sdss.sdss_flagval('ZWARNING', 'UNPLUGGED') | sdss.sdss_flagval('ZWARNING', 'BAD_TARGET')
            
            imain = np.where((wwhere(target,'*GALAXY*') | wwhere(target,'*GALAXY_BIG*') | wwhere(target,'*GALAXY_BRIGHT_CORE*')) & ((zans['ZWARNING_NOQSO'].data & bad_fiber) == 0))[0]
            nmain = len(imain)
            row['N_TARGET_MAIN'] = nmain
            if nmain > 0:
                row['SUCCESS_MAIN'] = 100.0 * len(np.where((zans[imain]['ZWARNING'].data == 0) & wwhere(zans[imain]['CLASS'].data, 'GALAXY*') | wwhere(zans[imain]['CLASS'].data, 'QSO*')))/ nmain

            ilrg1 = np.where((wwhere(target,'*GALAXY_RED*') | wwhere(target,'*GALAXY_BIG_II*') | wwhere(target, '*GAL_LOZ*') | wwhere(target, '*LRG')) & ((zans['ZWARNING_NOQSO'].data & bad_fiber) == 0))[0]
            nlrg1 = len(ilrg1)
            row['N_TARGET_LRG1'] = nlrg1
            if nlrg1 > 0:
                row['SUCCESS_LRG1'] = 100.0 * len(np.where((zans[ilrg1]['ZWARNING_NOQSO'].data == 0) & wwhere(zans[ilrg1]['CLASS'].data, 'GALAXY*')))/ nlrg1

            ilrg2 = np.where((wwhere(target,'*GAL_HIZ*') | wwhere(target,'*GAL_CMASS*')) & ((zans['ZWARNING_NOQSO'].data & bad_fiber) == 0))[0]
            nlrg2 = len(ilrg2)
            row['N_TARGET_LRG2'] = nlrg2
            if nlrg2 > 0:
                row['SUCCESS_LRG2'] = 100.0 * len(np.where((zans[ilrg2]['ZWARNING_NOQSO'].data == 0) & wwhere(zans[ilrg2]['CLASS'].data, 'GALAXY*')))/ nlrg2

            ielg = np.where((wwhere(target,'*ELG*')) & ((zans['ZWARNING_NOQSO'].data & bad_fiber) == 0))[0]
            nelg = len(ielg)
            row['N_TARGET_ELG'] = nelg
            if nelg > 0:
                row['SUCCESS_ELG'] = 100.0 * len(np.where((zans[ielg]['ZWARNING_NOQSO'].data == 0) & wwhere(zans[ielg]['CLASS'].data, 'GALAXY*')))/ nelg
            
            iqso = np.where((wwhere(target,'*QSO*')) & ((zans['ZWARNING'].data & bad_fiber) == 0))[0]
            nqso = len(iqso)
            row['N_TARGET_QSO'] = nqso
            if nqso > 0:
                row['SUCCESS_QSO'] = 100.0 * len(np.where((zans[iqso]['ZWARNING'].data == 0) & wwhere(zans[iqso]['CLASS'].data, 'QSO*')))/ nqso

            if row['STATUS1D'] == 'Done':
                return(row)

    if ptt.exists(spDiag1dlog):
        lastline = get_lastline(spDiag1dlog)
        if 'Successful completion' in lastline:
            #Case where this 1D log file completed, which is not a case that should ever occur
            row['STATUS1D'] = 'FAILED'
        else:
            #Case where this 1D log file isn't completed
            row['STATUS1D'] = 'RUNNING'
    else:
        row['STATUS1D'] = 'Pending'

    return(row)



def publicdata(row):
    global pulic_plate_data
    global spPlatelistMessage
    if pulic_plate_data is None:
        publicfile = ptt.join(getenv('SPECLOG_DIR'), 'opfiles', 'spPlateList.par')
        try:
            pulic_plate_data=read_table_yanny(publicfile,'SPPLATELIST')
            pulic_plate_data.convert_bytestring_to_unicode()
    
            #pulic_plate_data = Table(yanny(publicfile)['SPPLATELIST'])
        except:
            if spPlatelistMessage is False:
                splog.log('Empty of missing spPlateList.par file')
                spPlatelistMessage = True
            pulic_plate_data = None
            return(row)
    try:
        match = pulic_plate_data[np.where((pulic_plate_data['plate'].data == int(row['field'])) & (pulic_plate_data['mjd'].data == int(row['mjd'])))[0]]
    except:
        match = []
    if len(match) != 0:
        row['PUBLIC'] = match['public'][0]
        row['PLATEQUALITY'] = match['platequality'][0]
        row['QUALCOMMENTS'] = match['qualcomments'][0]
    return(row)


def get_survey(row, fieldfile):
    # get chunck info
    if int(row['FIELD']) < 16000:
        row = get_chunkinfo(row)
    else:
        if 'FIELD_CADENCE' not in row.keys():
            try:
                fmap = fits.getdata(fieldfile.replace('spField','spfibermap'),-1)
                #f"spfibermap-{row['FIELD']}-{row['MJD']}.fits", -1)
                row['FIELD_CADENCE'] = fmap[0]['FIELDCADENCE']
            except:
                pass
    
        if 'FIELD_CADENCE' in row.keys():
            if 'bright' in row['FIELD_CADENCE'].lower():
                row['SURVEY'] = 'mwm-bhm'
            elif 'dark' in row['FIELD_CADENCE'].lower():
                row['SURVEY'] = 'bhm-mwm'
            if 'dark_174x8' in row['FIELD_CADENCE'].lower() or 'dark_100x8' in row['FIELD_CADENCE'].lower():
                row['PROGRAMNAME'] = 'FPS-RM'
            else:
                row['PROGRAMNAME'] = 'FPS'
            
            if int(row['FIELD']) < 100000:
                if 'PROGRAMNAME' in row.keys():
                    if row['PROGRAMNAME'] == '':
                        row['PROGRAMNAME']  = 'FPS-COMMISSIONING'
                else:
                    row['PROGRAMNAME']  = 'FPS-COMMISSIONING'
                if 'SURVEY' in row.keys():
                    if row['SURVEY'] == '':
                        row['SURVEY']  = 'Commissioning'
                else:
                    row['SURVEY']  = 'Commissioning'
                if row['FIELD_CADENCE'] == '':
                    row['FIELD_CADENCE'] = 'Commissioning'


    return(row)

def get_2d_status(path,plan,row, epoch=False):

    yplan = yanny(ptt.join(path,plan))
    hdr = yplan.new_dict_from_pairs()
    row['OBSERVATORY'] = hdr['OBS']
    planlist = hdr['planfile2d']
    planlist = planlist.replace("'","").split(' ')
    logfile2d = []  # list of 2d log files that exists
    mjdlist = []
    statusdone = False
    statusrun  = False
    statusmissing = False
    st = []
    if epoch:
        path = ptt.join(path,'..')
        
    for plan2d in planlist:
        yplan_2d = yanny(ptt.join(path, plan2d))
        hdr_2d = yplan_2d.new_dict_from_pairs()
        mjdlist.append(hdr_2d['MJD'])
        thislogfile = plan2d.replace('spPlan2d', 'spDiag2d').replace('.par','.log')
        if ptt.exists(ptt.join(path,thislogfile)):
            logfile2d.append(thislogfile)
            lastline = get_lastline(ptt.join(path,thislogfile))
            if 'Successful completion' in lastline:
                #case where this 2d log file is completed
                statusdone = True
                st.append(1)
            else:
                #case where this wd log file is not completed
                statusrun = True
                st.append(0)
        else:
            #case where this 2d log file is missing
            statusmissing = True
            st.append(0)
    row['MJDLIST'] = ' '.join(mjdlist)
    if epoch:
        if sum(st) >0 and st[-1] ==1:
            statusmissing = False
            statusrun = False
        
    if statusmissing:
        row['STATUS2D'] = 'Pending'
    elif statusrun:
        row['STATUS2D'] = 'RUNNING'
    elif statusdone:
        row['STATUS2D'] = 'Done'
    
    return(row)

def get_lastline(filepath):
    with open(filepath, 'rb') as f:
        try:  # catch OSError in case of a one line file 
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()
    return(last_line)


def grep(filepath, grepstr):
    with open(filepath, 'rb', 0) as f:
        s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        if s.find(grepstr.encode('UTF-8')) != -1:
            return(True)
    return(False)



def get_cols(fieldfile, Field_list, run2d, run1d, legacy = False, skipcart=None, basehtml=None, epoch=False):
    if ptt.exists(fieldfile):
        thisrun1d = np.unique([ptt.basename(ptt.abspath(x)) for x in glob(ptt.join(ptt.dirname(fieldfile),'*')+'/')]).tolist()
        for dir_ in ['coadd','extraction','flat_extraction','epoch']:
            if dir_ in thisrun1d:
                thisrun1d.remove(dir_)
        for r1 in thisrun1d:
            row = {}
            if run1d is not None:
                if r1 not in run1d:
                    continue
            row['RUN1D'] = r1
        
    
            hdr = fits.getheader(fieldfile)

            if skipcart is not None:
                if hdr['CARTID'] in skipcart:
                    continue
            for col in Field_list.columns:
                try:
                    if type(hdr[col]) is str:
                        if len(hdr[col]) == 0: continue
                    row[col] = hdr[col]
                except:
                    continue
            row['RUN2D']         = run2d
            row['N_TOTAL']       = hdr['naxis2']
            try:
                row['FIELD_CADENCE'] = hdr['FIELDCAD']
            except:
                row['FIELD_CADENCE'] = ''
            try:
                row['SN2_G1']        = hdr['SPEC1_G']
                row['SN2_R1']        = hdr['SPEC1_R']
                row['SN2_I1']        = hdr['SPEC1_I']
            except:
                row['SN2_G1']        = np.nan
                row['SN2_R1']        = np.nan
                row['SN2_I1']        = np.nan
            try:
                row['SN2_G2']        = hdr['SPEC2_G']
                row['SN2_R2']        = hdr['SPEC2_R']
                row['SN2_I2']        = hdr['SPEC2_I']
            except:
                row['SN2_G2']        = np.nan
                row['SN2_R2']        = np.nan
                row['SN2_I2']        = np.nan
            try:
                row['MAPNAME']       = hdr['name']
            except:
                pass
            
            try:
                row['MOON_FRAC']    = hdr['MOONFRAC']
            except:
                row['MOON_FRAC']    = np.nan
            row['RACEN']         = hdr['RADEG']
            row['DECCEN']        = hdr['DECDEG']
            row['EPOCH']         = hdr['EQUINOX']
            row['STATUSCOMBINE'] = 'Done'

            if legacy:
                row['RACEN']         = hdr['RA']
                row['DECCEN']        = hdr['DEC']
                row['DERED_SN2_G1']  = hdr['SN2EXT1G']
                row['DERED_SN2_R1']  = hdr['SN2EXT1R']
                row['DERED_SN2_I1']  = hdr['SN2EXT1I']
                row['DERED_SN2_G2']  = hdr['SN2EXT2G']
                row['DERED_SN2_R2']  = hdr['SN2EXT2R']
                row['DERED_SN2_I2']  = hdr['SN2EXT2I']

            #get this from the file name since sometimes they are wrong in the file headers
            row['MJD']   = ptt.basename(fieldfile).replace('.fits','').split('-')[-1]
            row['FIELD'] = ptt.basename(fieldfile).replace('.fits','').split('-')[-2]

            row = get_survey(row, fieldfile)

            # Determine public data
            row = publicdata(row)

            row = getquality(row,basehtml, epoch=epoch)
            
            if epoch:
                combineplan = fieldfile.replace('spField', 'spPlancombepoch').replace('.fits','.par')
                row['STATUS2D'] = 'Done'
            else:
                combineplan = fieldfile.replace('spField', 'spPlancomb').replace('.fits','.par')
            row = get_2d_status(ptt.dirname(combineplan),ptt.basename(combineplan),row, epoch=epoch)
            if row['STATUS2D'] != 'RUNNING':
                row = read_spec1d(row, ptt.dirname(combineplan), fieldfile, epoch=epoch)
            elif row['STATUS2D'] != 'Done':
                row['STATUSCOMBINE'] = 'Pending'
                row['STATUS1D'] = 'Pending'
            if row['STATUS1D'] == 'Done':
                row = getoutputs(row,basehtml, epoch=epoch)
            else:
                row['PLOTS'] = ''
                row['DATA']  = ''
                
    
            for key in list(row.keys()):
                if key not in Field_list.columns:
                    row.pop(key)
            Field_list.add_row(row)
            
    else: ## no spField exists
        row={}
        row['RUN2D']         = run2d
        if epoch:
            thislogfile = fieldfile.replace('spField','spPlancombepoch').replace('.fits','.log')
        else:
            thislogfile = fieldfile.replace('spField','spDiagcomb').replace('.fits','.log')
        if ptt.exists(thislogfile):
            lastline = get_lastline(thislogfile)
            if 'Successful completion' in lastline:
                # case where spcombine completed but we are still missing spfield, 
                # so it must have failed
                row['STATUSCOMBINE'] = 'FAILED'
            else:
                # case where spcombine isn't completed
                abortline =  grep(thislogfile, 'ABORT')
                if abortline:
                    row['STATUSCOMBINE'] = 'FAILED'
                else:
                    row['STATUSCOMBINE'] = 'RUNNING'
        else:
            # case where spcombine log is missing
            row['STATUSCOMBINE'] = 'Pending'


        row['STATUS1D'] = 'Pending'
        row['DATA'] = ''
        row['PLOTS'] = ''
        row['SN2_G1']        = np.nan
        row['SN2_R1']        = np.nan
        row['SN2_I1']        = np.nan
        row['SN2_G2']        = np.nan
        row['SN2_R2']        = np.nan
        row['SN2_I2']        = np.nan
        row['FBADPIX']       = np.nan
        row['MOON_FRAC']     = np.nan
        row['EXPTIME']       = np.nan
        row['RACEN']         = np.nan
        row['DECCEN']        = np.nan
        row['FIELDSN2']      = np.nan
        row['FIELDQUALITY']  = 'bad'

        #get this from the file name since sometimes they are wrong in the file headers
        row['MJD']   = ptt.basename(fieldfile).replace('.fits','').split('-')[-1]
        row['FIELD'] = ptt.basename(fieldfile).replace('.fits','').split('-')[-2]
        row['PLOTSN'] = ''
        row = get_survey(row, fieldfile)
        # Determine public data
        row = publicdata(row)
        if epoch:
            combineplan = fieldfile.replace('spField', 'spPlancombepoch').replace('.fits','.par')
        else:
            combineplan = fieldfile.replace('spField', 'spPlancomb').replace('.fits','.par')
        
        row = get_2d_status(ptt.dirname(combineplan),ptt.basename(combineplan),row, epoch=epoch)

        
        if row['STATUSCOMBINE'] in ['RUNNING']:
            if row['STATUS2D'] == 'RUNNING':
                row['STATUS2D'] = 'FAILED'
                row['STATUSCOMBINE'] = 'Pending'
            else:
                thisrun1d = np.unique([ptt.basename(ptt.abspath(x)) for x in glob(ptt.join(ptt.dirname(fieldfile),'*')+'/')]).tolist()
                for dir_ in ['coadd','extraction','flat_extraction','epoch']:
                    if dir_ in thisrun1d:
                        thisrun1d.remove(dir_)
                for r1 in thisrun1d:
                    if run1d is not None:
                        if r1 not in run1d:
                            continue
                        spDiag1dlog = ptt.join(ptt.dirname(combineplan), r1, 'spDiag1d-'+row['FIELD']+'-'+row['MJD']+'.log')
                        if ptt.exists(spDiag1dlog):
                            row['STATUSCOMBINE'] = 'FAILED'
        elif row['STATUSCOMBINE'] == 'FAILED':
            if row['STATUS2D'] == 'RUNNING':
                row['STATUS2D'] = 'FAILED'
                row['STATUSCOMBINE'] = 'Pending'
                
        #row['RUN2D'] == run2d
            
        for key in list(row.keys()):
            if key not in Field_list.columns:
                row.pop(key)
        Field_list.add_row(row)
    return(Field_list)


def get_key(fp):
    filename = ptt.basename(fp)
    int_part = filename.split('-')[1]
    try:
        return int(int_part)
    except:
        return int(ptt.splitext(int_part)[0])

def fieldlist(create=False, topdir=getenv('BOSS_SPECTRO_REDUX'), run2d=[getenv('RUN2D')], run1d=[getenv('RUN1D')], outdir=None, 
              legacy=False, custom=None, skipcart=None, basehtml=None, datamodel= None, epoch=False, logfile=None, **kwrd):

    if datamodel is None:
        datamodel = ptt.join(getenv('IDLSPEC2D_DIR'), 'datamodel', 'fieldList_dm.par')

    global oplimit_filename
    oplimit_filename = ptt.join(getenv('IDLSPEC2D_DIR'),'examples','opLimits.par')

    if basehtml is None:
        basehtml = '../'
        if epoch is True:
            basehtml = '../../'

    if run1d is None: 
        run1d = run2d
    if skipcart is not None: 
        skipcart = np.atleast_1d(skipcart).astype(str).tolist()
    run1d = np.atleast_1d(run1d).astype(str).tolist()
    run2d = np.atleast_1d(run2d).astype(str).tolist()
    srun2d = '-'.join(run2d)
    if outdir is None:
        if epoch is True:
            outdir = ptt.join(topdir, srun2d, 'epoch')
        else:
            outdir = ptt.join(topdir, srun2d)
    makedirs(outdir, exist_ok = True)
    # if the create flag not set and the fieldlist file already exists then return the info in that file
    fitsfile = ptt.join(outdir, 'fieldlist-'+srun2d+'.fits')
    
    if logfile is None:
        splog.open(logfile = ptt.join(outdir, fitsfile.replace('.fits','.log')), backup=False)
        splog.log('Log file '+ptt.join(outdir, fitsfile.replace('.fits','.log'))+' opened '+ time.ctime())
    else:
        splog.open(logfile = logfile, backup=False)
        splog.log('Log file '+logfile+' opened '+ time.ctime())


    if ptt.exists(fitsfile) and create is False:
        return(Table(fits.getdata(fitsfile,1)))
    
    Field_list = Table(merge_dm(table=Table(), ext = 'FIELDLIST', name = 'FIELDLIST', dm =datamodel, splog=splog).data)
    Field_list.add_column(Column(name='PLOTSN', dtype=object))
    Field_list.add_column(Column(name='DATA', dtype=object))
    Field_list.add_column(Column(name='PLOTS', dtype=object))
    
    for r2 in run2d:
        path = ptt.join(topdir, r2)
        
        if epoch is not False:
            base = 'spPlancombepoch'
        elif custom is None:
            base = 'spPlancomb'
        if epoch is not False:
            fullfiles = sorted(glob(ptt.join(path,'*', 'epoch', base+'-*.par')), key=get_key)
        else:
            fullfiles = sorted(glob(ptt.join(path,'*', base+'-*.par')), key=get_key)

        for ff in fullfiles:
            splog.log('Reading '+ff)
            ff = ff.replace('.par','.fits').replace(base, 'spField')
            Field_list = get_cols(ff, Field_list, r2, run1d, legacy=legacy, skipcart=skipcart, basehtml=basehtml, epoch=epoch)



            
        if 'TILEID' in Field_list.colnames:
            #----------
            # Decide which fields constitute unique tiles with the required S/N,
            # then set QSURVEY=1.
            # Also insist that PROGNAME='main'.
            Field_list['QSURVEY'] = 0
            # First get the unique list of TILE
            tids = Field_list['TILEID'].data
            surv = np.char.lower(Field_list['SURVEY'].data.astype(str))
            fqual = np.char.lower(Field_list['FIELDQUALITY'].data.astype(str))
            tilelist = np.sort(np.unique(Field_list['TILEID'].data))
            qsurv = Field_list['QSURVEY']
            for itile in tilelist:
                indx = np.where((tids == itile) &
                                ((fqual == 'good') | (fqual == 'marginal')) &
                                ((surv  == 'bhm-mwm') | (surv  == 'bhm') | (surv  == 'mwm') | (surv  == 'boss')))[0]

                if (len(indx) > 0):
                    ibest = np.argmax(Field_list[indx]['FIELDSN2'].data)
                    qsurv[indx[ibest]] = 1


    if ptt.exists(ptt.join(outdir,'fieldlist-'+srun2d+'.fits')):
        remove(ptt.join(outdir,'fieldlist-'+srun2d+'.fits'))
    fhdu = merge_dm(table=Field_list, ext = 'FIELDLIST', name = 'FIELDLIST', dm =datamodel, splog=splog)
    plot_sky(outdir, Field_list, ptt.join(outdir,'fieldlist-'+srun2d+'.fits'))

    hdu = merge_dm(ext='Primary', hdr = {'RUN2D':srun2d,'Date':time.ctime()}, dm = datamodel, splog=splog)
                
    hdul = fits.HDUList([hdu, fhdu])
    hdul.writeto(ptt.join(outdir,'fieldlist-'+srun2d+'.fits'), overwrite=True)


    cols = {'FIELD':'FIELD','MJD':'MJD','OBS':'OBSERVATORY','PLOTS':'PLOTS','RACEN':'RACEN','DECCEN':'DECCEN',
            'RUN2D':'RUN2D','RUN1D':'RUN1D','DATA':'DATA','QUALITY':'FIELDQUALITY','EXPTIME':'EXPTIME',
            'SN^2':'FIELDSN2','N_gal':'N_GALAXY','N_QSO':'N_QSO','N_star':'N_STAR','N_unk':'N_UNKNOWN',
            'N_sky':'N_SKY','N_std':'N_STD','MOON_FRAC':'MOON_FRAC','SURVEY':'SURVEY','PROG':'PROGRAMNAME',
            'FIELD_CADENCE':'FIELD_CADENCE','DESIGNS':'DESIGNS','PUBLIC':'PUBLIC'}

    html_writer(basehtml, Field_list.copy(), outdir, 'fieldlist.html', srun2d, legacy, sort='field', order=cols, title = 'SDSS BOSS Spectroscopy Fields Observed List')
    html_writer(basehtml, Field_list.copy(), outdir, 'fieldlist-mjdsort.html', srun2d, legacy, sort='mjd',order=cols, title = 'SDSS BOSS Spectroscopy Fields Observed List')

    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'LCO'].copy(), outdir, 'fieldlist_LCO.html', srun2d, legacy, sort='field', order=cols, title = 'SDSS BOSS Spectroscopy LCO Fields Observed List')
    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'LCO'].copy(), outdir, 'fieldlist_LCO-mjdsort.html', srun2d, legacy, sort='mjd',order=cols, title = 'SDSS BOSS Spectroscopy LCO Fields Observed List')

    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'APO'].copy(), outdir, 'fieldlist_APO.html', srun2d, legacy, sort='field', order=cols, title = 'SDSS BOSS Spectroscopy APO Fields Observed List')
    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'APO'].copy(), outdir, 'fieldlist_APO-mjdsort.html', srun2d, legacy, sort='mjd',order=cols, title = 'SDSS BOSS Spectroscopy APO Fields Observed List')


    cols = {'FIELD':'FIELD','MJD':'MJD','OBS':'OBSERVATORY','PLOTS':'PLOTS','RACEN':'RACEN','DECCEN':'DECCEN',
            'RUN2D':'RUN2D','RUN1D':'RUN1D', 'SN2_G1':'SN2_G1','SN2_I1':'SN2_I1','SN2_G2':'SN2_G2',
            'SN2_I2':'SN2_I2','Badpix':'FBADPIX','SUCCESS_QSO':'SUCCESS_QSO','2D':'STATUS2D',
            'Combine':'STATUSCOMBINE','1D':'STATUS1D','SNPLOT':'PLOTSN','MOON_FRAC':'MOON_FRAC',
            'EXPTIME':'EXPTIME','QUALITY':'FIELDQUALITY','FIELD_CADENCE':'FIELD_CADENCE',
            'SURVEY':'SURVEY','PROG':'PROGRAMNAME','QUALCOMMENTS':'QUALCOMMENTS'}
            
    html_writer(basehtml, Field_list.copy(), outdir, 'fieldquality.html', srun2d, legacy, sort='field', order=cols, title='SDSS BOSS Spectroscopy Field Quality List')
    html_writer(basehtml, Field_list.copy(), outdir, 'fieldquality-mjdsort.html', srun2d, legacy, sort='mjd',order=cols, title='SDSS BOSS Spectroscopy Field Quality List')

    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'LCO'].copy(), outdir, 'fieldquality_LCO.html', srun2d, legacy, sort='field', order=cols, title='SDSS BOSS Spectroscopy LCO Field Quality List')
    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'LCO'].copy(), outdir, 'fieldquality_LCO-mjdsort.html', srun2d, legacy, sort='mjd',order=cols, title='SDSS BOSS Spectroscopy LCO Field Quality List')

    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'APO'].copy(), outdir, 'fieldquality_APO.html', srun2d, legacy, sort='field', order=cols, title='SDSS BOSS Spectroscopy APO Field Quality List')
    html_writer(basehtml, Field_list[Field_list['OBSERVATORY'] == 'APO'].copy(), outdir, 'fieldquality_APO-mjdsort.html', srun2d, legacy, sort='mjd',order=cols, title='SDSS BOSS Spectroscopy APO Field Quality List')


    splog.log('Successful completion of fieldlist at '+ time.ctime())
    splog.close()
    return(Field_list)

def color2hex(colorname):
    if colorname.strip().upper() == 'RED':    return('#FF0000')
    if colorname.strip().upper() == 'YELLOW': return('#909000')
    return('black')

def formatter(value, tl, strlimit):
    if strlimit is True:
        color = tl[tl['strval'].data == value]
    else:
        try:
            fval = float(value)
        except:
            return(value)
        color = tl[(tl['lovalue'].data <= fval) & (tl['hivalue'].data >= fval)]
    if len(color) == 1:
        color = color2hex(color['color'][0])
    else:
        color = 'black'
    return('<span style="color:'+color+';font-weight:bold;">'+ value + '</span>')

def html_format(column, cols_dic):
    global oplimits, oplimit_filename
    if oplimits is None:
        oplimits = yanny(oplimit_filename)

    raw_name = cols_dic[column.name]
    tl = Table(oplimits['TEXTLIMIT'])
    tl.convert_bytestring_to_unicode()
    match = tl[np.where((tl['field'].data == raw_name.upper()) & (tl['flavor'].data == 'SUMMARY') & (tl['camera'].data == '*'))[0]]
    strlimit = True
    if len(match) == 0:
        tl = Table(oplimits['SPECLIMIT'])
        tl.convert_bytestring_to_unicode()
        match = tl[np.where((tl['field'].data == raw_name.upper()) & (tl['flavor'].data == 'SUMMARY') & (tl['camera'].data == '*'))[0]]
        strlimit = False
    if len(match) != 0:
        column = column.apply(formatter, tl=match, strlimit=strlimit)
    return(column)

def html_writer(basehtml, Field_list, path, name, run2d, legacy, sort='field', order=None, title='SDSS Spectroscopy Fields Observed List'):
    run2d = run2d.replace('-',',')
    if order is not None:
        for col in order.keys():
            Field_list.rename_column(order[col],col)

        Field_list.keep_columns(list(order.keys()))
        Field_list = Field_list[list(order.keys())]
    fl_pd = Field_list.to_pandas()
    
    formats = {'RACEN':2,'DECCEN': 2,'SN2_G1':1,'SN2_I1':1,'SN2_G2':1,'SN2_I2': 1,'SN^2':1,
                 'Badpix':3,'SUCCESS_QSO':1,'MOON_FRAC':1,'EXPTIME':1}
    for key in formats.keys():
        if key in fl_pd.columns:
            if   formats[key] == 1:
                fl_pd[key] = fl_pd[key].map(lambda x: '%.1f' % float(x))
            elif formats[key] == 2:
                fl_pd[key] = fl_pd[key].map(lambda x: '%.2f' % float(x))
            elif formats[key] == 3:
                fl_pd[key] = fl_pd[key].map(lambda x: '%.3f' % float(x))
            if key in ['SN2_G1','SN2_I1','SN2_G2','SN2_I2','Badpix',
                       'MOON_FRAC','RACEN','DECCEN','EXPTIME','SN^2']:
                fl_pd[key] = fl_pd[key].map(lambda x: x.replace('nan',''))

    for i in fl_pd.index:
        if ('SN2_G1' in fl_pd.columns) & ('SN2_G2' in fl_pd.columns):
            sn1 = fl_pd.at[i, 'SN2_G1']
            sn2 = fl_pd.at[i, 'SN2_G2']
            fl_pd.at[i, 'SN2_G1'] = '' if sn2 != '0.0' else fl_pd.at[i, 'SN2_G1']
            fl_pd.at[i, 'SN2_G2'] = '' if sn1 != '0.0' else fl_pd.at[i, 'SN2_G2']
        if ('SN2_I1' in fl_pd.columns) & ('SN2_I2' in fl_pd.columns):
            sn1 = fl_pd.at[i, 'SN2_I1']
            sn2 = fl_pd.at[i, 'SN2_I2']
            fl_pd.at[i, 'SN2_I1'] = '' if sn2 != '0.0' else fl_pd.at[i, 'SN2_I1']
            fl_pd.at[i, 'SN2_I2'] = '' if sn1 != '0.0' else fl_pd.at[i, 'SN2_I2']
            if (fl_pd.at[i, 'SN2_I1'] == '0.0') & (fl_pd.at[i, 'SN2_I2'] == '0.0'):
                if fl_pd.at[i, 'SN2_G2'] == '': fl_pd.at[i, 'SN2_I2'] = ''
                if fl_pd.at[i, 'SN2_G1'] == '': fl_pd.at[i, 'SN2_I1'] = ''
            if (fl_pd.at[i, 'SN2_G1'] == '0.0') & (fl_pd.at[i, 'SN2_G2'] == '0.0'):
                if fl_pd.at[i, 'SN2_I2'] == '': fl_pd.at[i, 'SN2_G2'] = ''
                if fl_pd.at[i, 'SN2_I1'] == '': fl_pd.at[i, 'SN2_G1'] = ''
        if int(fl_pd.at[i,'FIELD']) >= 15000:
            clearcam = '1' if fl_pd.at[i,'OBS'] == 'LCO' else '2'
            for col in [f'SN2_I{clearcam}', f'SN2_G{clearcam}']:
                if col in fl_pd.columns:
                    fl_pd.loc[fl_pd.index[i],col] = fl_pd.loc[fl_pd.index[i],col].replace('0.0', '')

                    #fl_pd.at[i, c] = fl_pd.at[i, c].replace(['0', '0.0'], '')
#                if ('SN2_I1' in fl_pd.columns):
#                    fl_pd.at[i, 'SN2_I1'] = fl_pd.at[i, 'SN2_I1'].replace(['0', '0.0'], '')
#                if ('SN2_G1' in fl_pd.columns):
#                    fl_pd.at[i, 'SN2_G1'] = = fl_pd.at[i, 'SN2_I1'].replace(['0', '0.0'], '')
#            else:
#                if ('SN2_I2' in fl_pd.columns):
#                    fl_pd.at[i, 'SN2_I2'] = ''
#                if ('SN2_G2' in fl_pd.columns):
#                    fl_pd.at[i, 'SN2_G2'] = ''


    f2zero_cols = ['N_gal','N_QSO','N_star','N_unk','N_sky','N_std','SUCCESS_QSO']
    if '1D' in fl_pd.columns:
        idx = np.where(fl_pd['1D'].values != 'Done')
    else:
        idx = np.where(fl_pd['DATA'].values == '')
    for col in f2zero_cols:
        
        if col in fl_pd.columns:
            fl_pd[[col]] = fl_pd[[col]].astype(str)
            fl_pd.loc[fl_pd.index[idx],col] = fl_pd.loc[fl_pd.index[idx],col].replace('0', '')
            fl_pd.loc[fl_pd.index[idx],col] = fl_pd.loc[fl_pd.index[idx],col].replace('0.0', '')
            #fl_pd[[col]] = fl_pd[[col]].astype(str).replace(['0', '0.0'], '')

    f2s_cols = ['SN2_G1','SN2_I1','SN2_G2','SN2_I2', 'Badpix', 'SUCCESS_MAIN',
                'SUCCESS_LRG', 'SUCCESS_QSO', "%LRG1", "%LRG2", 'SN^2',
                'MOON_FRAC', 'EXPTIME','RACEN','DECCEN']
    for col in f2s_cols:
        if col in fl_pd.columns:
            fl_pd[[col]] = fl_pd[[col]].fillna('').astype(str)
    fl_pd = fl_pd.apply(html_format, cols_dic=order)


    
    if sort.lower() == 'mjd':
        fl_pd = fl_pd.sort_values(by=['MJD','FIELD'], ascending = False, key=lambda col: col.astype(int))
    else:
        fl_pd = fl_pd.sort_values(by=['FIELD','MJD'], key=lambda col: col.astype(int))
    html = fl_pd.to_html(escape=False, render_links=True, index=False)
    
    head = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
    <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
    <head>
        <title>"""+title+"""</title>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
        <style type="text/css">
            table { border-collapse: collapse; }
            th {
                padding: 2px;
                text-align: right;
                border: 1px solid black;
                font-weight: bold;
                }
            td {
                padding: 2px;
                text-align: right;
                border: 1px solid black;
                }
        </style>
     </head>"""
    mjd="{mjd:.3f}".format(mjd=float(astropy.time.Time(datetime.datetime.utcnow()).jd)-2400000.5)
    basehtml = 'href="'+basehtml
    if basehtml[-1] != '/':
        basehtml = basehtml+'/'
    basehtml = basehtml+run2d

    head2 = """
    <body>
    <h1>"""+title+"""</h1>
    <p>Last Update: """+time.ctime()+""", Last Update MJD: """+mjd+"""</p>
    <ul>
        <li><a """+basehtml+"""/">HOME</a></li>
        <li>Field list sorted by <a href="fieldlist.html">field</a>, <a href="fieldlist-mjdsort.html">MJD</a></li>
        <ul>
            <li>APO Field list sorted by <a href="fieldlist_APO.html">field</a>, <a href="fieldlist_APO-mjdsort.html">MJD</a></li>
            <li>LCO Field list sorted by <a href="fieldlist_LCO.html">field</a>, <a href="fieldlist_LCO-mjdsort.html">MJD</a></li>
        </ul>
        <li>Field quality sorted by <a href="fieldquality.html">field</a>, <a href="fieldquality-mjdsort.html">MJD</a></li>
        <ul>
            <li>APO Field quality sorted by <a href="fieldquality_APO.html">field</a>, <a href="fieldquality_APO-mjdsort.html">MJD</a></li>
            <li>LCO Field quality sorted by <a href="fieldquality_LCO.html">field</a>, <a href="fieldquality_LCO-mjdsort.html">MJD</a></li>
        </ul>
        <li>Field list as <a href="fieldlist-"""+run2d+""".fits">FITS</a></li>
    </ul>
    <p><a href="SDSSV.png"><img src="SDSSV_s.png"></a><a href="SDSSVc.png"><img src="SDSSVc_s.png"></a></p>
    <p><a href="SDSSV2.png"><img src="SDSSV2_s.png"></a><a href="SDSSV3.png"><img src="SDSSV3_s.png"></a></p>
    """
    if not legacy:
        head3="""<p>(S/N)^2 values are <b>not</b> corrected for galactic dust reddening</p>"""
    else:
        head3="""<p>(S/N)^2 values are corrected for galactic dust reddening</p>"""
    foot = """
    </body>
    </html>"""
    splog.log(ptt.join(path,name))

    with open(ptt.join(path,name), "w") as fhtml:
        fhtml.write(head)
        fhtml.write(head2)
        fhtml.write(head3)
        fhtml.write(html)
        fhtml.write(foot)
        fhtml.close()

if __name__ == '__main__' :
    """ 
    Build/load Fieldlist
    """
    parser = argparse.ArgumentParser(
            prog=ptt.basename(sys.argv[0]),
            description='Build/load BOSS Fieldlist')

    parser.add_argument('--create', '-c', action='store_true', help='Create Fieldlist')
    parser.add_argument('--topdir', type=str, help='Optional override value for the environment variable $BOSS_SPECTRO_REDUX', default = getenv('BOSS_SPECTRO_REDUX'))
    parser.add_argument('--run1d', type=str, help='Optional override value for the enviro variable $RUN1D', nargs='*', default=[getenv('RUN1D')])
    parser.add_argument('--run2d', type=str, help='Optional override value for the enviro variable $RUN2D', nargs='*', default=[getenv('RUN2D')])
    parser.add_argument('--outdir', type=str, help='Optional output directory (defaults to topdir/$RUN2D)', default=None) 
    parser.add_argument('--skipcart', type=str, help='Option list of cartridges to skip', nargs='*', default=None)
    parser.add_argument('--epoch', action='store_true', help='Produce FieldList for epoch coadds')
    parser.add_argument('--basehtml', type=str, help='html path for figure (defaults to relative from topdir)')
    parser.add_argument('--logfile', type=str, help='Manually Set logfile (including path)', default=None)
    parser.add_argument('--debug', action='store_true', help='Overrides the logger of the simplified error messages and prints standard python errors')
    args = parser.parse_args()

    splog.no_exception = args.debug
    Field_list = fieldlist(**vars(args))
