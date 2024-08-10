#!/usr/bin/env python3
from sdss_access.path import Path
from sdss_access import Access

from os import getenv
from astropy.table import Table
from pydl.pydlutils.yanny import read_table_yanny
import os.path as ptt
from glob import glob

class MissingFile(Exception):
    """Exception raise for missing plPlugMapM or confSummary file"""
    def __init__(self, file, filetype=None, message = None):
        if filetype is None:
            filetype = 'ConfSummary/ConfSummaryF'
        if message is None:
            message = f"{filetype} {ptt.basename(file)} is missing"
        self.message = message
        self.file = file
        super().__init__(message)
    


def find_plPlugMapM(mjd, plate, mapname, splog=None, release='sdsswork'):
    if getenv('SPECLOG_DIR') is None:
        if splog is None:
            print('ERROR: SPECLOG_DIR must be defined')
        else:
            splog.info('ERROR: SPECLOG_DIR must be defined')
        raise(MissingFile(None, message='SPECLOG_DIR must be defined as Environment variable'))
        exit()

    fibermap_file = None
    fibermap_file_list = []
    path   = Path(release=release, preserve_envvars=True)
    path_ops = {'mjd': mjd, 'plateid': int(plate)}
    try:
        fibermap_file = path.full('plPlugMapM', **path_ops)
        if path.exists('plPlugMapM', **path_ops):
            #local plPlugMapM file exists
            fibermap_file_list = [fibermap_file]
    except:
        fibermap_file = ptt.join(getenv('SPECLOG_DIR'),str(mjd),'plPlugMapM-'+mapname+'.par')
        fibermap_file_list = glob(fibermap_file)

    if len(fibermap_file_list) == 0:
        if splog is None:
            print('plPlugMapM file '+fibermap_file+' does not exists')
        else:
            splog.log('plPlugMapM file '+fibermap_file+' does not exists')
        raise(MissingFile(fibermap_file, filetype='plPlugMapM'))
        fibermap_file_list = [None]
    return(fibermap_file_list[0])


def find_confSummary(confid, obs=None, no_remote=False, splog=None, release='sdsswork'):
    path   = Path(release=release, preserve_envvars=True)
    if getenv('SDSSCORE_DIR') is None:
        if splog is None:
            print('ERROR: SDSSCORE_DIR must be defined')
        else:
            splog.info('ERROR: SDSSCORE_DIR must be defined')
        raise(MissingFile(None, message='SDSSCORE_DIR must be defined as Environment variable'))
        exit()
    if obs is None:
        obs = getenv('OBSERVATORY').lower()
    path_ops = {'configid': confid, 'obs': obs.lower()}
    fibermap_file = None

    if path.exists('confSummaryF_test', **path_ops):
        #local confSummaryF file exists
        fibermap_file = path.full('confSummaryF_test', **path_ops)
    elif path.exists('confSummary_test', **path_ops):
        #Local confSummary file exists
        fibermap_file = path.full('confSummary_test', **path_ops)
    else:
        #No confSummary file exists
        fibermap_file = path.full('confSummary', **path_ops)
        if splog is None:
            print('confSummary file '+fibermap_file+' does not exists')
        else:
            splog.info('confSummary file '+fibermap_file+' does not exists')
        raise(MissingFile(fibermap_file))
        fibermap_file = None
    return fibermap_file



def get_confSummary(confid, obs=None, no_remote=False, splog=None,
                    release='sdsswork', sort=False, filter=False):
    fibermap_file = find_confSummary(confid, obs=obs, no_remote=no_remote,
                                     splog=splog, release=release)

    if fibermap_file is None:
        return(Table())
    confSummary=read_table_yanny(fibermap_file,'FIBERMAP')
    if filter is True:
        confSummary=confSummary[np.where(confSummary['fiberType'] == 'BOSS')]
    if sort is True:
        confSummary.sort('fiberId')
    
    return(confSummary)
