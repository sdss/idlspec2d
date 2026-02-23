#!/usr/bin/env python3
from boss_drp.utils import putils
from boss_drp.utils.lock import lock, unlock
from boss_drp.utils.hash import create_hash
from boss_drp.sos.run_log2html import run_soslog2html

from pydl.pydlutils import yanny
from astropy.table import Table, unique
from astropy.io import fits
from os import getenv, remove, sep
from os import path as ptt
from glob import glob
import shutil
from collections import OrderedDict
import platform
from time import sleep
import re
import numpy as np
import os

try:
    import git
    if os.getenv('GIT_PYTHON_TRACE') is None:
        os.environ['GIT_PYTHON_TRACE'] = '2'  # 'full' gives even more details
except:
    git = None

def getLastMJD(silent=True):
    if ('sdss5' not in platform.node()) and (getenv('IDLSPEC2D_SOS', None) is None):
       print('mjd is required when not running at observatories')
       exit()
    else:
       path = ptt.join('/','data','spectro', '?????')
    def get_key(fp):
        if not ptt.isdir(fp): return(0)
        filename = ptt.basename(fp)
        int_part = filename.split()[0]
        return(int(int_part))
    files = sorted(glob(path),key=get_key)
    mjd = ptt.basename(files[-1])
    if silent is not True:
        print('Latest MJD %s' % mjd)
    return mjd

def read_sdHdrFix(sdHdrFix_file):
    try:
        return(yanny.read_table_yanny(sdHdrFix_file, 'OPHDRFIX'))
    except: 
        return(None)


def fixhdr(expid, hdrcards, mjd=None, obs=getenv('OBSERVATORY'), clobber=False, cameras='??', update=True, nogit=False):
    if mjd is None: mjd = getLastMJD()

    sdHdrFix_file = ptt.join(getenv('SDHDRFIX_DIR'), obs.lower(), 'sdHdrfix', 'sdHdrFix-'+str(mjd)+'.par')

    updates= None
    
    if clobber is not True: 
        updates = read_sdHdrFix(sdHdrFix_file)

    if updates is None:
        updates = Table(names=('fileroot', 'keyword', 'value'),
                        descriptions = ('Root of file name, without any ".fit" suffix', 'Keyword name', 'Keyword value (as a string)'),
                        dtype = ('S20', 'S9', 'S80'))

        updates.meta = OrderedDict({'MJD':  str(mjd)     +"   # Modified Julian Date for sdHdrFix file",
                                    'OBS':  str(obs)   +"     # Observatory" })
   
    for key in hdrcards.keys():
        frame = 'sdR-'+cameras+'-'+str(expid).zfill(8)
        updates.add_row((frame, key.upper(), hdrcards[key]))
        if (key.lower() == 'quality') & (update is True):
            fixSOSlog(frame,mjd,hdrcards[key],obs)
    updates = unique(updates, keys=['fileroot','keyword'], keep='last')

    print('Writing to: ',sdHdrFix_file)
    print(updates)
       
    if ptt.exists(sdHdrFix_file):
        remove(sdHdrFix_file)
    yanny.write_ndarray_to_yanny(sdHdrFix_file, updates, structnames='OPHDRFIX',hdr=updates.meta, comments=None)

    if (git is not None) and (not nogit):
        repo = git.Repo(getenv('SDHDRFIX_DIR'))  # Absolute path to repo
        relative_path = os.path.relpath(sdHdrFix_file, getenv('SDHDRFIX_DIR'))  # Convert to relative path
        repo.index.add([relative_path])
        print(f'Adding file to git repo')
    else:
        print('File not yet added to git repo... run the following commands to add it')
        print(f'cd {ptt.dirname(sdHdrFix_file)}')
        print(f'git add {ptt.basename(sdHdrFix_file)}')

    
def fixSOSlog(frame,mjd,quality,obs):
    logfiles = []
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sos',f'{mjd}',f'logfile-{mjd}.fits')))
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sosredo',f'{mjd}',f'logfile-{mjd}.fits')))
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sosredo','dev',f'{mjd}',f'logfile-{mjd}.fits')))
    for lf in logfiles:
        if ptt.exists(lf):
            if lock(f'{lf}', pause = 5):
                try:
                    with fits.open(lf, mode='update') as hdul:
                        print(f'Updating {lf}')
                        for ext in [1,2,3,4]:
                            try:
                                hdul[ext]
                            except:
                                continue
                            if '??' in frame:
                                ccds = ['b1','r1'] if obs.lower() == 'apo' else ['b2','r2']
                            else:
                                ccds = [None]

                            for ccd in ccds:
                                tframe = frame.replace('??',ccd) if ccd is not None else frame
                                if hdul[ext].data is None: continue
                                idx = np.where(hdul[ext].data['FILENAME'] == f'{tframe}.fit.gz')[0]
                                if len(idx) == 0:
                                    continue
                                else:
                                    hdul[ext].data['QUALITY'][idx[0]] = quality
                        hdul.flush()
                finally:
                    unlock(f'{lf}')
                    run_soslog2html(lf, mjd, obs)
                    
                test = create_hash(ptt.dirname(lf))
                if test:
                    print("\nsha1sum is locked")
            else:
                continue
        else:
            continue
class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __repr__(self):
        return '{{{0} - {1}}}'.format(self.start, self.end)

