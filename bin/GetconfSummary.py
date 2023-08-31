#!/usr/bin/env python3

from sdss_access.path import Path
from sdss_access import Access
from os import getenv, environ
from astropy.table import Table
from pydl.pydlutils.yanny import read_table_yanny, yanny, write_table_yanny
import os.path as ptt
from glob import glob


def find_plPlugMapM(mjd, plate, mapname, splog=None, release='sdsswork'):
    if getenv('SPECLOG_DIR') is None:
        splog.info('ERROR: SPECLOG_DIR must be defined')
        exit()
    fibermap_file = ptt.join(getenv('SPECLOG_DIR'),str(mjd),'plPlugMapM-'+mapname+'.par')
    fibermap_file_list = glob(fibermap_file)
    if len(fibermap_file_list) == 0:
        splog.info('plPlugMapM file '+fibermap_file+' does not exists')
        return(None)
    else:
        return(fibermap_file_list[0])
    
#
#    path   = Path(release='release', preserve_envvars=True)
#    print(path.lookup_names())
#    path_ops = {'mjd': mjd, 'plateid': int(plate)}
#    fibermap_file = path.full('plPlugMapM', **path_ops)
#    if path.exists('plPlugMapM', **path_ops):
#        pass
#    elif path.exists('plPlugMapM', **path_ops, remote=True):
#        access = Access(release='release')
#        access.remote()
#        access.add('plPlugMapM', **path_ops)
#        access.set_stream()
#        valid = access.commit()
#        if valid is False:
#            splog.log('plPlugMapM file '+fibermap_file+' does not exists')
#            return(None)
#    return(fibermap_file)


def find_confSummary(confid, obs=None, no_remote=False, splog=None, release='sdsswork'):
    path   = Path(release=release, preserve_envvars=True)

    if obs is None:
        obs = getenv('OBSERVATORY').lower()
    path_ops = {'configid': confid, 'obs': obs.lower()}
    fibermap_file = None
    #print(path.full('confSummaryF_test', **path_ops))
    if path.exists('confSummaryF_test', **path_ops):
        #local confSummaryF file exists
        fibermap_file = path.full('confSummaryF_test', **path_ops)
    elif (not no_remote):
        access = Access(release=release)
        if (path.exists('confSummaryF_test', **path_ops, remote=True)):
            #Remote confSummaryF file exists
            fibermap_file = path.full('confSummaryF_test', **path_ops)
            access.remote()
            access.add('confSummaryF_test', **path_ops)
            access.set_stream()
            valid = access.commit()
            if valid is False:
                fibermap_file = None
    if fibermap_file is None:
        if path.exists('confSummary_test', **path_ops):
            #Local confSummary file exists
            fibermap_file = path.full('confSummary_test', **path_ops)
        elif (not no_remote):
            access = Access(release=release)
            if (path.exists('confSummary_test', **path_ops, remote=True)):
                #Remote confSummary file exists
                fibermap_file = path.full('confSummary_test', **path_ops)
                access.remote()
                access.add('confSummary_test', **path_ops)
                access.set_stream()
                valid = access.commit()
                if valid is False:
                    if splog is None:
                        print('confSummary file '+fibermap_file+' does not exists')
                    else:
                        splog.info('confSummary file '+fibermap_file+' does not exists')
                    return None
    if fibermap_file is None:
        #No confSummary file exists
        fibermap_file = path.full('confSummary', **path_ops)
        if splog is None:
            print('confSummary file '+fibermap_file+' does not exists')
        else:
            splog.info('confSummary file '+fibermap_file+' does not exists')
        return None
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
