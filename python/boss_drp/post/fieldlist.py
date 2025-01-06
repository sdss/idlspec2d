#!/usr/bin/env python3
from boss_drp import idlspec2d_dir, favicon
from boss_drp.field import field_to_string, Field
from boss_drp.utils import match as wwhere
from boss_drp.utils import (grep, get_lastline, merge_dm, Splog, jdate, retry)
from boss_drp.post import plot_sky_targets, plot_sky_locations

import argparse
import sys
import os
import os.path as ptt
import numpy as np
from pydl.pydlutils.yanny import yanny, read_table_yanny
from pydl.pydlutils import sdss
from astropy.io import fits
from astropy.table import Table, Column, unique
import astropy.time
from glob import glob
import time
import datetime
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import gc
from jinja2 import Template

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
    sdss.set_maskbits(maskbits_file=ptt.join(os.getenv("IDLUTILS_DIR"),"data","sdss","sdssMaskbits.par"))
except:
    splog.log('Environmental Varable IDLUTILS_DIR must be set')
    exit()

def getquality(row, dereddened_sn2=False, rawsn2=False):
    field = Field(fieldlist_name.basehtml, row['RUN2D'], row['FIELD'],
                  epoch=fieldlist_name.epoch)
    sfield = field.field_str
    plotsn = ptt.join(field.dir(), 'spSN2d-'+sfield+'-'+str(row['MJD'])+'.pdf')
    if not ptt.exists(plotsn):
        plotsn = plotsn.replace('.pdf','.ps')
    plotsn = ptt.join(ptt.relpath(field.dir(), fieldlist_name.name), ptt.basename(plotsn))
    plotsn = path_to_html(plotsn)
    row['PLOTSN'] = '<a href="'+plotsn+'">SNPLOT</a>'
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
        row['FIELDSN2'] = np.nanmin(fieldsn2[valid])
        
        if dereddened_sn2:
            deredsn2 = np.array([DERED_SN2_G1,DERED_SN2_I1,DERED_SN2_G2,DERED_SN2_I2])
            row['DEREDSN2'] = np.nanmin(deredsn2[valid])
        
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

def getoutputs(row, field):
#    field = Field(fieldlist_name.basehtml, row['RUN2D'], row['FIELD'],
#                  custom_name=fieldlist_name.custom_name,
#                  epoch = fieldlist_name.epoch)
    PLOTS = field.png_dir(row['RUN1D'],row['MJD'], pathbase = fieldlist_name.basehtml)
    PLOTS = ptt.join(ptt.relpath(ptt.dirname(PLOTS), fieldlist_name.name), ptt.basename(PLOTS))
    PLOTS = path_to_html(PLOTS, dir=True)

    DATA  = field.spec_dir(row['MJD'], pathbase = fieldlist_name.basehtml)
    DATA = ptt.join(ptt.relpath(ptt.dirname(DATA), fieldlist_name.name), ptt.basename(DATA))
    DATA = path_to_html(DATA, dir=True)
    
    row['PLOTS'] = '<a href="'+PLOTS+'">PLOTS</a>'
    row['DATA'] = '<a href="'+ DATA+ '">DATA</a>'
    return(row)
    
def get_chunkinfo(row):
    global chunkdata
    if chunkdata is None:
        chunkfile = ptt.join(os.getenv('PLATELIST_DIR'), 'platePlans.par')
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



def read_spec1d(row, field_class):#path, fieldfile):
    spZfile     = ptt.join(field_class.dir(), row['RUN1D'],
                            'spZbest-'+row['FIELD']+'-'+row['MJD']+'.fits')
    spDiag1dlog = ptt.join(field_class.dir(), row['RUN1D'],
                            'spDiag1d-'+row['FIELD']+'-'+row['MJD']+'.log')
    if ptt.exists(spZfile):
        strmjd = str(row['MJD']).strip()
        run2d = row['RUN2D']
        run1d = row['RUN1D']
        
        try:
            zans = Table(fits.getdata(spZfile, 1))
            plug = Table(fits.getdata(field_class.spField,5))
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
        publicfile = ptt.join(os.getenv('SPECLOG_DIR'), 'opfiles', 'spPlateList.par')
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
                fmap = None
                pass
        else:
            fmap = None
    
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

    if 'DESIGN_VERS' not in row.keys():
        try:
            if fmap is None:
                fmap = fits.getdata(fieldfile.replace('spField','spfibermap'),-1)
            row['DESIGN_VERS'] = fmap[0]['DESIGN_VERS']
        except:
            pass

    return(row)


def get_DesignMode(field_class, row): #path,plan,row):
    designMode = []

    path = field_class.dir()
    if fieldlist_name.epoch:
        path = ptt.join(path,'..')
        
    for plan2d in field_class.plan2d:
        yp2d = yanny(ptt.join(path,plan2d))
        hdr = yp2d.new_dict_from_pairs()
        mjd = hdr['MJD']
        field = hdr['fieldname']
        spFibermap = 'spfibermap-'+field_to_string(field)+'-'+str(mjd)+'.fits'
        if ptt.exists(ptt.join(path,spFibermap)):
            with fits.open(ptt.join(path,spFibermap)) as hdul:
                for i, ext in enumerate(hdul[1].data['EXTNAME']):
                    en = (ext.replace('.par','').replace('plPlugMapM-','')
                                               .replace('confSummaryF-','')
                                               .replace('confSummary-',''))
                    if en in field_class.configIDs:
                        try:
                            designMode.append(hdul[1].data['DESIGN_MODE'][i])
                        except Exception as e:
                            print(e)
                            pass
    if len(designMode) == 0:
        row['DESIGN_MODE'] = ''
    else:
        designMode = np.unique(np.asarray(designMode))
        if len(designMode) > 1:
            row['DESIGN_MODE'] = 'mixed_mode'
        else:
            row['DESIGN_MODE'] = designMode[0]
    return(row)
        

def get_2d_status(field_class, row):#path,plan,row):

    yplan = yanny(ptt.join(field_class.dir(),field_class.plancomb))
    hdr = yplan.new_dict_from_pairs()
    try:
        row['OBSERVATORY'] = hdr['OBS']
    except:
        if yplan['SPEXP']['name'][0][0].astype(str).split('-')[1] in ['b2','r2']:
            row['OBSERVATORY'] = 'LCO'
        else:
            row['OBSERVATORY'] = 'APO'
    planlist = hdr['planfile2d']
    planlist = planlist.replace("'","").split(' ')
    logfile2d = []  # list of 2d log files that exists
    mjdlist = []
    statusdone = False
    statusrun  = False
    statusmissing = False
    st = []
    if fieldlist_name.epoch:
        path = ptt.join(field_class.dir(),'..')
    else:
        path = field_class.dir()
        
    field_class.plan2d = planlist
    field_class.configIDs = yplan['SPEXP']['mapname'].astype(str)
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
    if fieldlist_name.epoch:
        if sum(st) >0 and (st[-1] ==1 or int(mjdlist[-1]) < jdate.mjd - 2):
            statusmissing = False
            statusrun = False
        
    if statusmissing:
        row['STATUS2D'] = 'Pending'
    elif statusrun:
        row['STATUS2D'] = 'RUNNING'
    elif statusdone:
        row['STATUS2D'] = 'Done'
    
    return(row)


def get_cols(field_class, Field_list, run2d, run1d, legacy = False, skipcart=None):
    if ptt.exists(field_class.spField):
        thisrun1d = np.unique([ptt.basename(ptt.abspath(x)) for x in glob(ptt.join(ptt.dirname(field_class.spField),'*')+'/')]).tolist()
        for dir_ in ['coadd','extraction','flat_extraction','epoch']:
            if dir_ in thisrun1d:
                thisrun1d.remove(dir_)
        for r1 in thisrun1d:
            row = {}
            if run1d is not None:
                if r1 not in run1d:
                    continue
            row['RUN1D'] = r1
        
    
            hdr = fits.getheader(field_class.spField)

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
            row['MJD']   = ptt.basename(field_class.spField).replace('.fits','').split('-')[-1]
            row['FIELD'] = ptt.basename(field_class.spField).replace('.fits','').split('-')[-2]
        
            field_class.run2d = row['RUN2D']
            field_class.run1d = row['RUN1D']
            field_class.field = row['FIELD']
            field_class.mjd   = row['MJD']
            field_class.set()
            
            
            
            row = get_survey(row, field_class.spField)

            # Determine public data
            row = publicdata(row)
            row = getquality(row)
            
            if fieldlist_name.epoch:
                row['STATUS2D'] = 'Done'
            row = get_2d_status(field_class, row)
            row = get_DesignMode(field_class, row)
            if row['STATUS2D'] != 'RUNNING':
                row = read_spec1d(row, field_class)
            elif row['STATUS2D'] != 'Done':
                row['STATUSCOMBINE'] = 'Pending'
                row['STATUS1D'] = 'Pending'
            if row['STATUS1D'] == 'Done':
                row = getoutputs(row, field_class)
            else:
                row['PLOTS'] = ''
                row['DATA']  = ''
                
    
            for key in list(row.keys()):
                if key not in Field_list.columns:
                    row.pop(key)
            Field_list.add_row(row)
            
    else: ## no spField exists
        splog.info(f'{ptt.basename(field_class.spField)} not found, checking intermediate status')
        row={}
        row['RUN2D']         = run2d
        if fieldlist_name.epoch:
            thislogfile = field_class.spField.replace('spField','spPlancombepoch').replace('.fits','.log')
        else:
            thislogfile = field_class.spField.replace('spField','spDiagcomb').replace('.fits','.log')
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
        row['MJD']   = ptt.basename(field_class.spField).replace('.fits','').split('-')[-1]
        row['FIELD'] = ptt.basename(field_class.spField).replace('.fits','').split('-')[-2]
        row['PLOTSN'] = ''
        
        field_class.run2d = row['RUN2D']
        field_class.run1d = None
        field_class.field = row['FIELD']
        field_class.mjd   = row['MJD']
        field_class.set()
        row = get_survey(row, field_class.spField)
        # Determine public data
        row = publicdata(row)
        row = get_2d_status(field_class, row)
        row = get_DesignMode(field_class, row)

        
        if row['STATUSCOMBINE'] in ['RUNNING']:
            if row['STATUS2D'] == 'RUNNING':
                row['STATUS2D'] = 'FAILED'
                row['STATUSCOMBINE'] = 'Pending'
            else:
                thisrun1d = np.unique([ptt.basename(ptt.abspath(x)) for x in glob(ptt.join(ptt.dirname(field_class.spField),'*')+'/')]).tolist()
                for dir_ in ['coadd','extraction','flat_extraction','epoch']:
                    if dir_ in thisrun1d:
                        thisrun1d.remove(dir_)
                for r1 in thisrun1d:
                    if run1d is not None:
                        if r1 not in run1d:
                            continue
                        spDiag1dlog = ptt.join(field_class.dir(), r1, 'spDiag1d-'+row['FIELD']+'-'+row['MJD']+'.log')
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

def fieldlist(create=False, topdir=os.getenv('BOSS_SPECTRO_REDUX'), run2d=[os.getenv('RUN2D')],
              run1d=[os.getenv('RUN1D')], outdir=None, legacy=False, custom=None,
              skipcart=None, basehtml=None, datamodel= None, epoch=False, return_tab=False,
              logfile=None, noplot=False, field=None, mjd=None, debug=False, **kwrd):
              
    if datamodel is None:
        datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'fieldList_dm.par')

    global oplimit_filename
    oplimit_filename = ptt.join(idlspec2d_dir,'examples','opLimits.par')

    if basehtml is None:
        basehtml = topdir if topdir is not None else  os.getenv('BOSS_SPECTRO_REDUX')
    fieldlist_name.basehtml = basehtml
#    if (field is not None) and (mjd is not None):
#        if basehtml is None:
#            basehtml = '../'
#            if epoch is True:
#                basehtml = '../../'
#    else:
#        basehtml = '../../../'

    if run1d is None: 
        run1d = run2d
    if skipcart is not None: 
        skipcart = np.atleast_1d(skipcart).astype(str).tolist()
    run1d = np.atleast_1d(run1d).astype(str).tolist()
    run2d = np.atleast_1d(run2d).astype(str).tolist()
    srun2d = '-'.join(run2d)
    
    fieldlist_name.build(topdir, srun2d, epoch=epoch, custom_name = custom,
                             logfile=logfile, outdir=outdir)
    os.makedirs(fieldlist_name.outdir, exist_ok = True)
    # if the create flag not set and the fieldlist file already exists then return the info in that file
    fitsfile = fieldlist_name.name
    if (field is None) and (mjd is None):
        global splog
        if logfile is None:
            tmpext = '.tmp'
            splog.open(logfile = ptt.join(outdir, fitsfile.replace('.fits','.log')), backup=False)
            splog.log('Log file '+ptt.join(outdir, fitsfile.replace('.fits','.log'))+' opened '+ time.ctime())
        else:
            tmpext = ptt.basename(logfile).replace('.log','.tmp').replace('fieldlist-','')
            splog.open(logfile = logfile, backup=False)
            splog.log('Log file '+logfile+' opened '+ time.ctime())
    else:
        splog = kwrd['fmsplog']

    splog.no_exception = debug
    if ptt.exists(fitsfile) and create is False:
        return(Table(fits.getdata(fitsfile,1)))
    
    Field_list = Table(merge_dm(table=Table(), ext = 'FIELDLIST', name = 'FIELDLIST', dm =datamodel, splog=splog).data)
    Field_list.add_column(Column(name='PLOTSN', dtype=object))
    Field_list.add_column(Column(name='DATA', dtype=object))
    Field_list.add_column(Column(name='PLOTS', dtype=object))
    
    if (field is not None) and (mjd is not None):
        if ptt.exists(fitsfile):
            Field_list = Table(fits.getdata(fitsfile))
            idx  = np.where((Field_list['FIELD'] == field) & (Field_list['MJD'] == int(mjd)))[0]
            if len(idx) > 0:
                Field_list.remove_rows(idx)
 
 
    for r2 in run2d:
        path = ptt.join(topdir, r2)
        
        if epoch is not False:
            base = 'spPlancombepoch'
        elif custom is None:
            base = 'spPlancomb'
        field_class = Field(topdir, r2, '*', epoch=epoch)
        fullfiles = sorted(glob(ptt.join(field_class.dir(), base+'-*.par')), key=get_key)

        nfields = len(fullfiles)
        for ifield, ff in enumerate(fullfiles):
            if (field is not None) and (mjd is not None):
                if f"{field}-{mjd}" not in ff:
                    continue
            field_class.plancomb = ff
            ff = ff.replace('.par','.fits').replace(base, 'spField')
            field_class.spField = ff
            splog.log('Reading '+ff+ f' ({ifield+1}/{nfields})')
            Field_list = get_cols(field_class, Field_list, r2, run1d,legacy=legacy,
                                  skipcart=skipcart)



            
        if 'TILEID' in Field_list.colnames:
            splog.log('Checking for best field in each unique tile')
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
            ibest = None
            indx = None
            for itile in tilelist:
                indx = np.where((tids == itile) &
                                ((fqual == 'good') | (fqual == 'marginal')) &
                                ((surv  == 'bhm-mwm') | (surv  == 'bhm') | (surv  == 'mwm') | (surv  == 'boss')))[0]

                if (len(indx) > 0):
                    ibest = np.argmax(Field_list[indx]['FIELDSN2'].data)
                    qsurv[indx[ibest]] = 1
            del tids
            del surv
            del fqual
            del tilelist
            del qsurv
            del indx
            try:
                del ibest
            except:
                pass
        del fullfiles

    if (field is None) and (mjd is None):
        write_fieldlist(Field_list, srun2d, datamodel, legacy=legacy, noplot=noplot)

    if (field is not None) and (mjd is not None):
        idx  = np.where((Field_list['FIELD'] == field) & (Field_list['MJD'] == int(mjd)))[0]
        if len(idx) > 0:
            return(Field_list[idx[0]])
 

    if not noplot:
        if not return_tab:
            del Field_list
            Field_list = None
        else:
            Field_list = None
        retry(plot_sky_locations, retries=3, delay=5, logger=splog.log)
        if not ptt.exists(ptt.join(fieldlist_name.outdir,'SDSSV2.png')):
            summary_names.set(topdir, srun2d,epoch=epoch, custom=custom)
            retry(plot_sky_targets, retries=3, delay=5, logger=splog.log, nobs=True)
    elif return_tab:
        Field_list = Table.read(fitsfile)
        Field_list.convert_bytestring_to_unicode()
    else:
        Field_list = None
    splog.log('Successful completion of fieldlist at '+ time.ctime())
    splog.close()
    return(Field_list)

def write_fieldlist(Field_list, srun2d, datamodel, legacy=False, noplot=False):
    cols = {'FIELD':'FIELD','MJD':'MJD','OBSERVATORY':'OBS','PLOTS':'PLOTS','RACEN':'RACEN','DECCEN':'DECCEN',
            'RUN2D':'RUN2D','RUN1D':'RUN1D','DATA':'DATA','FIELDQUALITY':'QUALITY','EXPTIME':'EXPTIME',
            'FIELDSN2':'SN^2','N_GALAXY':'N_gal','N_QSO':'N_QSO','N_STAR':'N_star','N_UNKNOWN':'N_unk',
            'N_SKY':'N_sky','N_STD':'N_std','MOON_FRAC':'MOON_FRAC','SURVEY':'SURVEY','PROGRAMNAME':'PROG',
            'DESIGN_MODE':'DESIGN_MODE',#'DESIGN_VERS':'DESIGN_VERS',
            'FIELD_CADENCE':'FIELD_CADENCE','DESIGNS':'DESIGNS','PUBLIC':'PUBLIC'}
    try:
        html_writer(Field_list, 'fieldlist', srun2d, legacy,
                    order=cols, title = 'SDSS BOSS Spectroscopy {obs} Fields Observed List')
    except:
        time.sleep(60)
        html_writer(Field_list,'fieldlist', srun2d, legacy,
                    order=cols, title = 'SDSS BOSS Spectroscopy {obs} Fields Observed List')
    cols = {'FIELD':'FIELD','MJD':'MJD','OBSERVATORY':'OBS','PLOTS':'PLOTS','RACEN':'RACEN','DECCEN':'DECCEN',
            'RUN2D':'RUN2D','RUN1D':'RUN1D', 'SN2_G1':'SN2_G1','SN2_I1':'SN2_I1','SN2_G2':'SN2_G2',
            'SN2_I2':'SN2_I2','FBADPIX':'Badpix','SUCCESS_QSO':'SUCCESS_QSO','STATUS2D':'2D',
            'STATUSCOMBINE':'Combine','STATUS1D':'1D','PLOTSN':'SNPLOT','MOON_FRAC':'MOON_FRAC',
            'EXPTIME':'EXPTIME','FIELDQUALITY':'QUALITY','FIELD_CADENCE':'FIELD_CADENCE',
            'DESIGN_MODE':'DESIGN_MODE',#'DESIGN_VERS':'DESIGN_VERS',
            'SURVEY':'SURVEY','PROGRAMNAME':'PROG','QUALCOMMENTS':'QUALCOMMENTS'}

    try:
        html_writer(Field_list, 'fieldquality', srun2d, legacy,
                    order=cols, title='SDSS BOSS Spectroscopy {obs} Field Quality List')
    except:
        time.sleep(60)
        html_writer(Field_list,'fieldquality', srun2d, legacy,
                    order=cols, title='SDSS BOSS Spectroscopy {obs} Field Quality List')

    splog.info('Formatting Fits File')
    Field_list = merge_dm(table=Field_list, ext = 'FIELDLIST', name = 'FIELDLIST',
                    dm =datamodel, splog=splog, drop_cols=['PLOTSN','DATA','PLOTS'])

    hdu = merge_dm(ext='Primary', hdr = {'RUN2D':srun2d,'Date':time.ctime()}, dm = datamodel, splog=splog)
       
    splog.info('writing: '+fieldlist_name.name)
    hdul = fits.HDUList([hdu, Field_list])
    try:
        hdul.writeto(fieldlist_name.name+fieldlist_name.tmpext, overwrite=True)
        os.rename(fieldlist_name.name+fieldlist_name.tmpext, fieldlist_name.name)
    except:
        time.sleep(60)
        hdul.writeto(fieldlist_name.name+fieldlist_name.tmpext, overwrite=True)
        os.rename(fieldlist_name.name+fieldlist_name.tmpext, fieldlist_name.name)

    Field_list = None
    hdul = None
    hdu = None
    return

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

    raw_name = column.name #cols_dic[column.name]
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

def html_writer(Field_list, name, run2d, legacy, sorts=['field','mjd'], order=None,
                title='SDSS Spectroscopy Fields Observed List', obss=[None, 'LCO','APO']):
    #run2d = run2d.replace('-',',')
    if order is not None:
        #for col in order.keys():
            #Field_list.rename_column(order[col],col)

        #Field_list.keep_columns(list(order.keys()))
        #Field_list = Field_list[list(order.keys())]
        #fl_pd = Field_list[list(order.keys())].to_pandas()
        fl_pd = Field_list[list(order.keys())].to_pandas()
        fl_pd = fl_pd.rename(columns=order)
    else:
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

    mjd="{mjd:.3f}".format(mjd=jdate.mjd)
        
    basehtml = ptt.relpath(ptt.join(fieldlist_name.basehtml, run2d),
                           ptt.join(fieldlist_name.outdir, fieldlist_name.name))
    basehtml = path_to_html(fieldlist_name.basehtml, dir=True)

    #basehtml = basehtml
    #if basehtml[-1] != '/':
    #    basehtml = basehtml+'/'
    #basehtml = basehtml+run2d

    red = '<b>not </b>' if not legacy else ''
    foot = """
    </body>
    </html>"""

    fl_pd_full = fl_pd.copy()
    fl_pd = None
    for sort in sorts:
        if sort.lower() == 'mjd':
            fl_pd_full = fl_pd_full.sort_values(by=['MJD','FIELD'], ascending = False, key=lambda col: col.astype(int))
        else:
            fl_pd_full = fl_pd_full.sort_values(by=['FIELD','MJD'], key=lambda col: col.astype(int))
        for obs in obss:
            fl_pd = fl_pd_full.copy()
            
            tname = fieldlist_name.html[name+'_mjd'] if sort.lower() == 'mjd' else fieldlist_name.html[name]
            

            if obs is not None:
                fl_pd = fl_pd.loc[fl_pd['OBS'] == obs]
                obsstr = obs.upper()
                fos = '_'+obs.upper()
            else:
                obsstr = ''
                fos = ''
            tname = tname.format(obs=fos)
            
            fl_pd.columns = fl_pd.columns.str.replace('_',' ', regex=False)
            html = fl_pd.to_html(escape=False, render_links=True, index=False)
    

            splog.log(ptt.join(fieldlist_name.outdir,tname))
            template = ptt.join(idlspec2d_dir,'templates','html','fieldlist_template.html')
            
            jinja_data = dict(RUN2D=run2d,date=time.ctime(),MJD=mjd,obs=obsstr,
                            favicon=favicon, basehtml=basehtml,red=red,
                            FIELDLIST_TABLE=html, title=title.format(obs=obsstr))
            with open(ptt.join(fieldlist_name.outdir,tname), "w", encoding="utf-8") as output_file:
                with open(template) as template_file:
                    j2_template = Template(template_file.read())
                    output_file.write(j2_template.render(jinja_data))
                

                        

    fl_pd = None
    fl_pd_full = None
    head = head2 = head3 = html = foot = thead2 = None


