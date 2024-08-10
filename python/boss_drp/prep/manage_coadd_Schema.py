#!/usr/bin/env python3
from boss_drp.utils import Splog

import numpy as np
from os import getenv, makedirs
import os.path as ptt
import sys
from astropy.table import Table, vstack
import astropy.time
import argparse
from collections import OrderedDict
from pydl.pydlutils import yanny
from subprocess import getoutput
from datetime import date
import pandas as pd
from time import ctime

length_checker = np.vectorize(len)


def clean_tab(Tab):
    for col in Tab.colnames: 
        if Tab[col].dtype == 'bool':
            Tab.rename_column(col, col+'_BOOL')
            new_col = np.zeros(len(Tab), dtype=int)
            new_col[np.where(Tab[col+'_BOOL'].data)[0]] =1
            Tab[col] = new_col
            Tab.remove_column(col+'_BOOL')
        elif Tab[col].dtype == 'O':
            Tab.rename_column(col, col+'_OBJ')
            Tab[col] = [[""]]
            Tab.remove_column(col+'_OBJ')
        elif Tab[col].dtype == int: pass
        elif Tab[col].dtype == float: pass
        elif type(Tab[col].data) == np.ndarray:
            #print(col, str(len(max(Tab[col].data.ravel()))))
            #print(Tab[col].dtype)
            #print(length_checker(Tab[col].data).max())
            Tab[col]=Tab[col].astype('|S'+str(length_checker(Tab[col].data).max()))
            #Tab[col]=Tab[col].astype('|S'+str(len(max(Tab[col].data.ravel()))))
    return(Tab)

def restore_tab(Tab):
    for col in Tab.colnames:
        if Tab[col].dtype == int: 
            Tab.rename_column(col, col+'_BOOL')
            Tab[col] = False
            new_col = np.full(len(Tab), False, dtype=bool)
            new_col[np.where(Tab[col+'_BOOL'].data==1)[0]] = True
            Tab[col]=new_col
            Tab.remove_column(col+'_BOOL')
        else:
            Tab[col].astype(str)
    return(Tab)


def matchd(arr1, arr2):
    if len(arr2.shape) == 1: arr2 = np.atleast_2d(arr2)
    app = arr2.shape[1] - arr1.shape[1]
    if app < 0: app = 0
    arr1 = np.pad(arr1, [(0,0),(0, app)], mode='constant', constant_values='')
    app = arr1.shape[1] - arr2.shape[1]
    if app < 0: app = 0
    arr2 = np.pad(arr2, [(0,0),(0,app)], mode='constant', constant_values='')
    return(arr1,np.atleast_1d(arr2.squeeze()))


def where2d1d(arr2d, arr1d):
    merged = pd.concat([pd.DataFrame(arr2d),pd.DataFrame([arr1d],index=['c'])])
    return(np.array(merged.duplicated(keep='last').tolist()[:-1]))
    
    
def find_row(Tab, row):
    tcart, rcart = matchd(Tab['CARTON'].data.astype(str), row['CARTON'])
    tcat,  rcat  = matchd(Tab['SDSS_ID'].data.astype(str), row['SDSS_ID'])
    tmjd,  rmjd  = matchd(Tab['MJD'].data.astype(str), row['MJD'])
    tprog, rprog = matchd(Tab['PROGRAM'].data.astype(str), row['PROGRAM'])
    tleg, rleg = matchd(Tab['LEGACY'].data.astype(str), row['LEGACY'])
    
    i = np.where((Tab['NAME'].data.astype(str) == row['NAME']) & 
                 (where2d1d(tcart,rcart)) &
                 (where2d1d(tcat,rcat)) &
                 (where2d1d(tmjd,rmjd)) &
                 (where2d1d(tprog,rprog)) &
                 (where2d1d(tleg,rleg)) &
                 (Tab['CADENCE'].data  == row['CADENCE']) &
                 (Tab['DR'].data == row['DR']) &
                 (Tab['RERUN1D'].data == row['RERUN1D']))[0]
    return(i)


def merge(Tab1, Tab2):

    for col in Tab1.colnames:
        if len(Tab1[col].shape) > 1:
            c1 = Tab1[col].data
            c2 = Tab2[col].data
            app = c2.shape[1] - c1.shape[1]
            if app < 0: app = 0  
            c1 = np.pad(c1, [(0,0),(0, app)],mode='constant', constant_values='')
            app = c1.shape[1] - c2.shape[1]
            if app < 0: app = 0
            c2 = np.pad(c2, [(0,0),(0, app)],mode='constant', constant_values='')

            Tab1[col] = c1
            Tab2[col] = c2
        elif Tab1[col].dtype == bool: continue
        elif np.issubdtype(Tab1[col].dtype, np.number): continue
        else:
            c1 = Tab1[col].data
            c2 = Tab2[col].data
            slen1 = length_checker(c1).max()
#int(str(c1.dtype).replace('|S',''))
            slen2 =  length_checker(c1).max() #int(str(c2.dtype).replace('|S',''))
            slen = max(slen1,slen2)
            c1 = c1.astype('|S'+str(slen))
            c2 = c2.astype('|S'+str(slen))
            Tab1[col] = c1
            Tab2[col] = c2
    return(vstack([Tab1,Tab2]))

def manage_coadd_Schema(Name, topdir=None, run2d=None, DR=False, CARTON=None, CATID=None, RERUN1D=False, 
                        CADENCE=0.0, MJD=None, PROGRAM=None, ACTIVE=True, legacy=None, show=False, 
                        coaddfile=None, survey='BHM', SDSS_GEN='SDSS-V', use_catid=False,
                        use_firstcarton=False):


    if CARTON is None: CARTON = [""]
    if CATID is None: CATID = [""]
    if MJD is None: MJD = [""]
    if PROGRAM is None: PROGRAM = [""]
    if legacy is None: legacy = [""]
    if topdir is None: topdir = getenv('BOSS_SPECTRO_REDUX')
    if run2d is None: run2d = getenv('RUN2D')

    if coaddfile is None:
        coaddfile = ptt.join(topdir,run2d,'fields','SDSSV_BHM_COADDS.par')

    if show is True:
        if ptt.exists(coaddfile):
            print(restore_tab(yanny.read_table_yanny(coaddfile, 'SCHEMA')))
        else:
            print(f'{coaddfile} is missing')
        return


    caf = ptt.join(topdir, run2d, coaddfile)
    print(caf)
    if (DR is False) & (CARTON == [""]) & (CATID == [""]) & (RERUN1D is False) & (CADENCE == 0) & (MJD == [""]):
        if ptt.exists(caf):
            COADDS = restore_tab(yanny.read_table_yanny(caf, 'SCHEMA'))
    else:
        newCoadd = Table({'NAME':[Name],'DR':[DR],'CARTON':[CARTON],'SDSS_ID':[CATID], 'LEGACY': [legacy],
                          'RERUN1D':[RERUN1D], 'CADENCE':[CADENCE], 'MJD':[MJD], 'ACTIVE':[ACTIVE],
                          'USE_CATID':[use_catid], 'USE_FIRSTCARTON':[use_firstcarton],
                          'PROGRAM':[PROGRAM]})

        if ptt.exists(caf):
             COADDS = restore_tab(yanny.read_table_yanny(caf, 'SCHEMA'))
        else: COADDS = newCoadd

        if  ptt.exists(caf) :
            i = find_row(COADDS, newCoadd[0])
            if len(i) != 0 :
                COADDS.remove_rows(i)
                COADDS = merge(COADDS, newCoadd)
            else:
                i = find_row(COADDS, newCoadd)
                if len(i) == 0: COADDS = merge(COADDS, newCoadd)
        else: COADDS = newCoadd
    print(COADDS)
    COADDS = clean_tab(COADDS)
    COADDS.meta=OrderedDict({'#':'',
                             '# Schema Description for custom BOSS Coadds':'',
                             '# ':'',
                             '# Last Updated '+ctime():'',
                             '# Written by manage_coadd_Schema':'',
                             '#  ': '',
                             '#   ': 'Legacy is a currently unutilized, but there for future versions to include tags from SDSS-IV,-III etc in addition to current tag',
                             '#    ': '',
                             'SDSS_GEN': SDSS_GEN  + '  # Generation of SDSS',
                             'SURVEY':   survey    + '  # Survey/Mapper',
                             'Filename': coaddfile + '  # COADD schema filename (this file)',
                             })
    COADDS.convert_unicode_to_bytestring()
    yanny.write_table_yanny(COADDS, caf, tablename = 'SCHEMA', overwrite=True)

