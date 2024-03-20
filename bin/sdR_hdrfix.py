#!/usr/bin/env python

from pathlib import Path
from pydl.pydlutils import yanny
from astropy.table import Table, vstack
from os import getenv, environ, remove, sep, symlink, unlink
from os import path as ptt
from sys import argv
import argparse
from astropy.io import fits
from astropy.table import Table, Column, unique
from glob import glob
from collections import OrderedDict
import numpy as np
from multiprocessing import Process
import subprocess
from astropy.time import Time
from time import sleep
import platform
from termcolor import colored
from run_soslog2html import run_soslog2html

def getLastMJD(silent=True):
    if 'sdss5' not in platform.node():
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


def fixhdr(expid, hdrcards, mjd=None, obs=getenv('OBSERVATORY'), clobber=False, cameras='??'):
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
        if key.lower() == 'quality':
            fixSOSlog(frame,mjd,hdrcards[key],obs)
    updates = unique(updates, keys=['fileroot','keyword'], keep='last')

    print('Writing to: ',sdHdrFix_file)
    print(updates)
       
    if ptt.exists(sdHdrFix_file):
        remove(sdHdrFix_file)
    yanny.write_ndarray_to_yanny(sdHdrFix_file, updates, structnames='OPHDRFIX',hdr=updates.meta, comments=None)

def fixSOSlog(frame,mjd,quality,obs):
    logfiles = []
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sos',f'{mjd}',f'logfile-{mjd}.fits')))
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sosredo',f'{mjd}',f'logfile-{mjd}.fits')))
    logfiles.append(ptt.abspath(ptt.join(sep,'data','boss','sosredo','dev',f'{mjd}',f'logfile-{mjd}.fits')))
    for lf in logfiles:
        if ptt.exists(lf):
            try:
                symlink(f'{lf}',f'{lf}.lock')
            except:
                while ptt.exists(f'{lf}.lock'):
                    sleep(5)
                symlink(f'{lf}',f'{lf}.lock')
        else:
            continue
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
                        
                    idx = np.where(hdul[ext].data['FILENAME'] == f'{tframe}.fit.gz')[0]
                    if len(idx) == 0:
                        continue
                    else:
                        hdul[ext].data['QUALITY'][idx[0]] = quality
            hdul.flush()
        unlink(f'{lf}.lock')
        run_soslog2html(lf, mjd, obs)


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __repr__(self):
        return '{{{0} - {1}}}'.format(self.start, self.end)

if __name__ == '__main__' :
    msg = "\n    "+colored("At current only use if still exposing or don't run SOS after (skip and note in Night Log (and/or email) if uncertain)",'red')


    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=60)

    parser = argparse.ArgumentParser(
            prog=ptt.basename(argv[0]),
            description='',
            formatter_class=formatter,
            epilog='one or more update options are required')
    
    parser.add_argument('expid', help='Exposure ID')
    parser.add_argument('--mjd', '-m', help = 'mjd of file (default: '+str(getLastMJD(silent=True))+')', required=False, default = getLastMJD(silent=True))
    if 'sdss5' not in platform.node():
        if getenv('OBSERVATORY') is None: 
            obs_req = True
            req = parser.add_argument_group('Required arguments')
            req.add_argument('--obs', help = 'Set Observatory', required=True, choices=['APO','LCO'])
        else: 
            obs_req = False
            parser.add_argument('--obs', help = 'Set Observatory (default:'+getenv('OBSERVATORY')+')', default=getenv('OBSERVATORY'), required=False, choices=['APO','LCO'])

    parser.add_argument('--clobber', help='clobber sdHdrFix file', action='store_true')

    if 'sdss5' in platform.node():
        cam = ['b1','r1','??'] if getenv('OBSERVATORY').lower == 'apo' else ['b2','r2','??']
        carts = ['FPS-N'] if getenv('OBSERVATORY').lower == 'apo' else ['FPS-S']
    else:
        cam = ['b1','b2','r1','r2','??']
        carts = ['FPS-S', 'FPS-N']

    parser.add_argument('--cameras', type=str.lower, choices=cam, default='??', help='Cameras for hdr update (?? for all cameras) [default:??]')

    qualgroup = parser.add_argument_group("Optional Quality Update (exclusive)"+msg)
    qualgroup.add_argument('--bad', '-b', help = 'flag as quality=bad', action = 'store_true')
    qualgroup.add_argument('--test', '-t', help = 'flag as quality=test', action = 'store_true')

    lampgroup = parser.add_argument_group('Optional lamp/screen keys to Update (1:on, 0:off)')
    lampgroup.add_argument('--FF',      type=str, nargs=4, choices=['0','1'], help='Flat Field Lamp')
    lampgroup.add_argument('--FFS',     type=str, nargs=8, choices=['0','1'], help='Flat Field Screen')
    lampgroup.add_argument('--NE',      type=str, nargs=4, choices=['0','1'], help='Ne arc lamp')
    if 'sdss5' in platform.node():
        if getenv('OBSERVATORY').lower == 'apo':
            lampgroup.add_argument('--HGCD',    type=str, nargs=4, choices=['0','1'], help='HgCd arc Lamp')
        else:
            lampgroup.add_argument('--HEAR',    type=str, nargs=4, choices=['0','1'], help='HeAr arc Lamp')
    else:
        lampgroup.add_argument('--HGCD',    type=str, nargs=4, choices=['0','1'], help='HeCd arc Lamp')
        lampgroup.add_argument('--HEAR',    type=str, nargs=4, choices=['0','1'], help='HeAr arc Lamp')
    lampgroup.add_argument('--arc',       help = 'short cut to set all relevant arc lamps to 1 1 1 1', action='store_true')
    lampgroup.add_argument('--flat',      help = 'short cut to set FF = 1 1 1 1 & FFS =  1 1 1 1 1 1 1 1 ', action = 'store_true')
    lampgroup.add_argument('--hartmann',type=str.capitalize, choices = ['Out', 'Right', 'Left', 'Closed'], help='Hartmann Door Status')

    commongroup = parser.add_argument_group("Optional Common keys to Update"+msg)
    commongroup.add_argument('--quality', type=str.lower, choices=['excellent', 'test', 'bad'], help = 'Set Quality flat of exposures')

    specialgroup = parser.add_argument_group("Optional Specialized Keys to Update \n    "+colored(" Only use if still exposing or don't run SOS after (skip and note in Night Log (and/or email) if uncertain)",'red'))
    specialgroup.add_argument('--flavor',  type=str.lower, choices=['bias','dark','flat','arc','science','smear'], help='Type/Flavor of exposure')
    specialgroup.add_argument('--exptime', type=float, help='Exposure length (s)')
    specialgroup.add_argument('--tai-beg', type=float, help='Starting time (tai) of exposure')
    specialgroup.add_argument('--cartid',  type=str.upper, choices=carts, help='Cartridge Mounted')
    specialgroup.add_argument('--fieldid', type=int, metavar='FIELDID',  help="FieldID")
    specialgroup.add_argument('--confid',  type=int, metavar='CONFIGID', help="ConfigureID")
    specialgroup.add_argument('--designid',type=int, metavar='DESIGNID', help='DesignID')

    keygroup = parser.add_argument_group("Manually update a key \n    "+colored(" Only use if still exposing or don't run SOS after (skip and note in Night Log (and/or email) if uncertain)",'red'))
    keygroup.add_argument('--key', '-k', help = 'header keyword to update (required if value is set)', default=None)
    keygroup.add_argument('--value', '-v', help = 'updated header keyword value (required if key is set)', default=None)

    args =  parser.parse_args()
    

    if 'sdss5' in platform.node():
        args.obs = getenv('OBSERVATORY')


    if args.arc:
        if args.obs.lower == 'apo':
            args.NE = ['1','1','1','1']
            args.HGCD = ['1','1','1','1']
        else:
            args.NE = ['1','1','1','1']
            args.HEAR = ['1','1','1','1']
    if args.flat:
            args.FF  = ['1','1','1','1']
            args.FFS = ['1','1','1','1','1','1','1','1']
 
    if args.key is not None:
        if args.key.upper in ['EXPOSURE', 'MJD']:
            print('ERROR: Invalid Keys. Exiting')
            exit()    

    updates = {} 
   
    if args.bad: 
        updates = {'quality': 'bad'}
    elif args.test:
        updates = {'quality': 'test'}
    else:
        for key in vars(args).keys():
            if key in ['expid', 'mjd', 'obs', 'clobber', 'bad', 'test', 'key', 
                       'value','rerun', 'cameras', 'observer', 'arc', 'flat']: continue
            if vars(args)[key] is None: continue
            if type(vars(args)[key]) == list:
                updates[key] = ' '.join(vars(args)[key])
            else: 
                updates[key] = str(vars(args)[key])
        if args.key is not None and args.value is not None: 
            updates[args.key] = args.value

    
    for key in list(updates.keys()):
        if key in ['fieldid', 'confid', 'designid', 'flavor', 'exptime', 'tai-beg', 'cartid']:
            confirm = ''
            while confirm.lower() != 'y':
                confirm = input('Do you really want to edit '+key+'? (y/[n]) ')
                if confirm.lower() == '': confirm = 'n'
                if confirm.lower() == 'n': 
                    updates.pop(key)
                    print('Skipping '+key)
                    break
                
    if len(updates) == 0: 
        parser.print_help()    
    else:
        fixhdr(args.expid, updates, mjd = args.mjd, obs = args.obs, clobber=args.clobber, 
               cameras=args.cameras)

