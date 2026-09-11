#!/usr/bin/env python3
from boss_drp import MOUNTAIN
from boss_drp.utils.lock import lock, unlock
from boss_drp.utils.hash import create_hash
from boss_drp.sos.run_log2html import run_soslog2html
from boss_drp.utils.splog import splog

from pydl.pydlutils import yanny
from astropy.table import Table, unique
from astropy.io import fits
from os import getenv, remove, sep
from os import path as ptt
from glob import glob
from collections import OrderedDict
import numpy as np
import os
from pathlib import Path
import re

try:
    import git
    if os.getenv('GIT_PYTHON_TRACE') is None:
        os.environ['GIT_PYTHON_TRACE'] = '2'  # 'full' gives even more details
except:
    git = None

def getLastMJD(silent=True, print_func=print):
    if not MOUNTAIN:
       print_func('mjd is required when not running at observatories')
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
        print_func('Latest MJD %s' % mjd)
    return mjd

def read_sdHdrFix(sdHdrFix_file):
    try:
        return(yanny.read_table_yanny(sdHdrFix_file, 'OPHDRFIX'))
    except: 
        return(None)


def fixhdr(expid, hdrcards, mjd=None, obs=getenv('OBSERVATORY'), clobber=False, 
           cameras='??', update=True, nogit=False, print_func = print):
    if mjd is None: mjd = getLastMJD(print_func=print_func)

    sdHdrFix_file = ptt.join(getenv('SDHDRFIX_DIR'), obs.lower(), 'sdHdrfix', 'sdHdrFix-'+str(mjd)+'.par')
    if lock(sdHdrFix_file, pause=5, niter=12):
        try:
            updates= None
            
            if clobber is not True: 
                updates = read_sdHdrFix(sdHdrFix_file)

            if updates is None:
                updates = Table(names=('fileroot', 'keyword', 'value'),
                                descriptions = ('Root of file name, without any ".fit" suffix', 
                                                'Keyword name', 'Keyword value (as a string)'),
                                dtype = ('S20', 'S9', 'S80'))

                updates.meta = OrderedDict({'MJD': str(mjd) +"   # Modified Julian Date for sdHdrFix file",
                                            'OBS': str(obs) +"     # Observatory" })
        
            for key in hdrcards.keys():
                frame = 'sdR-'+cameras+'-'+str(expid).zfill(8)
                updates.add_row((frame, key.upper(), hdrcards[key]))
                if (key.lower() == 'quality') & (update is True):
                    fixSOSlog(frame,mjd,hdrcards[key],obs, print_func=print_func)
            updates = unique(updates, keys=['fileroot','keyword'], keep='last')

            print_func('Writing to: ',sdHdrFix_file)
            print_func(updates)
            
            if ptt.exists(sdHdrFix_file):
                remove(sdHdrFix_file)
            yanny.write_ndarray_to_yanny(sdHdrFix_file, updates, structnames='OPHDRFIX',
                                         hdr=updates.meta, comments=None)
        finally:
            unlock(sdHdrFix_file)

    if (git is not None) and (not nogit):
        repo = git.Repo(getenv('SDHDRFIX_DIR'))  # Absolute path to repo
        relative_path = os.path.relpath(sdHdrFix_file, getenv('SDHDRFIX_DIR'))  # Convert to relative path
        repo.index.add([relative_path])
        print_func(f'Adding file to git repo')
    else:
        print_func('File not yet added to git repo... run the following commands to add it')
        print_func(f'cd {ptt.dirname(sdHdrFix_file)}')
        print_func(f'git add {ptt.basename(sdHdrFix_file)}')

    
def fixSOSlog(frame,mjd,quality,obs, print_func=print):
    logfiles = []
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sos',f'{mjd}',f'logfile-{mjd}.fits')))
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sosredo',f'{mjd}',f'logfile-{mjd}.fits')))
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sosredo','dev',f'{mjd}',f'logfile-{mjd}.fits')))
    for lf in logfiles:
        if ptt.exists(lf):
            if lock(f'{lf}', pause = 5):
                try:
                    with fits.open(lf, mode='update') as hdul:
                        print_runc(f'Updating {lf}')
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
                    print_func("\nsha1sum is locked")
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



_failures = {'SCIENCE':['ABORT: Reject science as too bright'],
             'FLAT':[r'WARNING: Reject flat as too faint',
                     r'WARNING: Reject flat: \d+(?:\.\d+)?% bad pixels',
                     r'WARNING: Reject flat: \d+(?:\.\d+)? saturated rows'],
             'ARC':[r'WARNING: Reject arc: \d+(?:\.\d+)?% bad pixels',
                    r'WARNING: Reject arc: \d+(?:\.\d+)? saturated rows']
}

_aborts = {'SCIENCE':r'ABORT: Unable to reduce science exposure',
           'FLAT':r'ABORT: Unable to reduce flat',
           'ARC':r'ABORT: Unable to reduce arc'}

_ignore = {'SCIENCE':[r'ABORT: Reject science: Flat-field screens are closed!',
                      r'ABORT: Reject science: Flat-field lamps turned on!',
                      r'WARNING: Hartmann doors closed'],
          'FLAT':[r'WARNING: Reject flat: Flat-field lamps not turned on!',
                  r'WARNING: Reject flat: Flat-field screens not closed!',
                  r'WARNING: Hartmann doors closed'],
          'ARC':[r'WARNING: \d+/\d+ (?:Ne|HgCd|HeAr) lamps are off',
                 r'WARNING: Reject arc: Neither Ne nor HeAr lamps turned on!',
                 r'WARNING: Reject arc: Neither Ne nor HgCd lamps are on!',
                 r'WARNING: Reject Arc: Flat-field lamps turned on!',
                 r'WARNING: Hartmann doors closed',
                 r'WARNING: Reject arc: Flat-field screens not closed!']}
  
def flag_bad(mjd, sosdir='/data/boss/sos/', exposure=None,
             obs=getenv('OBSERVATORY', 'APO'), flavor=None, 
             print_only=False):

    logfile = Path(sosdir) / f'{mjd}' / f'logfile-{mjd}.fits'

    if MOUNTAIN:
        if lock(logfile, pause=5, niter=12):
            try:
                if not logfile.exists():
                    return
                with fits.open(logfile) as hdul:
                    message = hdul[5].data
                if message is None: 
                    return
            finally:
                unlock(logfile)
    else:
        print_only = True
        if not logfile.exists():
            return
        with fits.open(logfile) as hdul:
            message = hdul[5].data
        if message is None: 
            return    

    if exposure is not None:
        message = message[message['EXPNUM'] == exposure]

    if flavor is None:
        flavors = ['SCIENCE', 'ARC', 'FLAT']
    else:
        flavors = [flavor]

    fix_hdr_args = dict(clobber=False, update= False, nogit=True, print_func = splog.info)
    for expnum in set(message['EXPNUM']):
        _mess = message[message['EXPNUM'] == expnum]
        # Determine which flavor, if any, has an abort.
        abort_flavors = []

        for _flavor in flavors:
            for text in _mess['TEXT']:
                if isinstance(text, bytes):
                    text = text.decode()

                if re.search(_aborts[_flavor], text):
                    abort_flavors.append(_flavor)
                    break

        # If there is an abort, use that flavor.
        if abort_flavors:
            _flavor = abort_flavors[0]
        else:
            continue
        if obs== 'LCO':
            cameras = set(['b2','r2'])
        elif obs == 'APO' and int(mjd) < 59145:
            cameras = set(['b1','r1','b2','r2'])
        else:
            cameras = set(['b1','r1'])
        # cameras = set(_mess['CAMERA'])
        failed_cameras = set()

        # Check failure messages for the selected flavors.
        for _cam in cameras:
            mess = _mess[_mess['CAMERA'] == _cam]

            # Check whether this camera should be ignored
            ignore = False
            for row in mess:
                text = row['TEXT']
                if isinstance(text, bytes):
                    text = text.decode()

                if any(
                    re.search(pattern, text)
                    for _flavor in flavors
                    for pattern in _ignore[_flavor]
                ):
                    ignore = True
                    break

            if ignore:
                continue

            for row in mess:
                text = row['TEXT']
                if isinstance(text, bytes):
                    text = text.decode()

                matches = [
                    pattern
                    for _flavor in flavors
                    for pattern in _failures[_flavor]
                    if re.search(pattern, text)
                ]

                if matches:
                    failed_cameras.add(_cam)
                    break
        if not failed_cameras:
            continue
        # All cameras failed -> use ??
        if failed_cameras == cameras:
            if not print_only:
                fixhdr( expnum, {'quality': 'bad'}, obs, camera='??', **fix_hdr_args)
            splog.info(f'OPHDRFIX sdR-??-{int(expnum):08d} QUALITY bad')

        # Only some cameras failed -> flag individually
        else:
            for _cam in failed_cameras:
                if not print_only:
                    fixhdr( expnum, {'quality': 'bad'}, obs, camera=_cam, **fix_hdr_args)
                splog.info(f'OPHDRFIX sdR-{_cam}-{int(expnum):08d} QUALITY bad')

if __name__ == "__main__":
    from glob import glob
    for mjd_ in glob(str(Path(getenv('BOSS_SOS_S'))/'?????')): 
        mjd_ = Path(mjd_).name
        flag_bad(mjd_,getenv('BOSS_SOS_S'),obs='LCO',print_only=True)

    for mjd_ in glob(str(Path(getenv('BOSS_SOS_N'))/'?????')): 
        mjd_ = Path(mjd_).name
        flag_bad(mjd_,getenv('BOSS_SOS_N'),obs='APO',print_only=True)



#LCO 60099