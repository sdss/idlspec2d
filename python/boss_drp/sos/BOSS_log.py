#!/usr/bin/env python3
from boss_drp.utils import jdate

import os
import os.path as ptt
import sys
import argparse
import pandas as pd
from astropy.io import fits
from glob import glob
import astropy.time
import numpy as np
import platform
try:
    from termcolor import colored
except:
    def colored(text,color):
        print(text)
import time
import datetime
from pydl.pydlutils import yanny

import builtins
try:
    from ansi2html import Ansi2HTMLConverter
    from bs4 import BeautifulSoup
except:
    pass
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib
import re

class type_ref_point:
    def __init__(self, obs, mjd, expid, exptime, blue, red, WAVEMID_blue = 0, WAVEMID_red = 0):
        self.obs = obs
        self.mjd = mjd
        self.expid = expid
        self.exptime = exptime
        self.f_ref = {}
        self.WAVEMID_ref = {}
        if self.obs.lower() == 'apo':
            self.f_ref['b1'] = blue
            self.f_ref['r1'] = red
            self.WAVEMID_ref['b1'] = WAVEMID_blue
            self.WAVEMID_ref['r1'] = WAVEMID_red
        elif self.obs.lower() == 'lco':
            self.f_ref['b2'] = blue
            self.f_ref['r2'] = red
            self.WAVEMID_ref['b2'] = WAVEMID_blue
            self.WAVEMID_ref['r2'] = WAVEMID_red
class ref_point:
    def __init__(self, obs, flat, arc, dark, bias):
        self.obs  = obs
        self.flat = flat
        self.arc  = arc
        self.dark = dark
        self.bias = bias


output = []

flat = type_ref_point('lco', 60002, '00009216', 100, 800.485595703125,29655.01953125) #median
arc  = type_ref_point('lco', 60002, '00009217', 45 , 6.0288896560668945,20.390592575073242,
                      WAVEMID_blue = 4995.053639483913, WAVEMID_red = 8019.644330955543) #median
bias = type_ref_point('lco', 60002, '00009210', 0,   4.296058654785156,5.612983226776123) #98%
dark = type_ref_point('lco', 60002, '00009210', 300, 4.6704792976379395,5.790493965148925) #98%
lco_ref = ref_point('lco', flat, arc, dark, bias)
                           

flat = type_ref_point('apo', 60000, '00353056', 25.1,6001.1640625,31622.7578125)
arc  = type_ref_point('apo', 60000, '00353032', 4.09,8.71207046508789,29.698558807373047,
                      WAVEMID_blue = 4903.772623256195, WAVEMID_red = 8009.896220483651)
bias = type_ref_point('apo', 60000, '00353021', 0,3.7101967334747314,5.478268146514893)
dark = type_ref_point('apo', 60000, '00353022', 20,3.681579351425171,5.54632043838501)
apo_ref = ref_point('apo', flat, arc, dark, bias)

def get_key(fp):
    filename = ptt.splitext(ptt.splitext(ptt.basename(fp))[0])[0]
    int_part = filename.split('-')[2]
    try:
        return int(int_part)
    except:
        return int(ptt.splitext(int_part)[0])

def get_SOS_log_val(SOS_log, ffile, ext, col):
    if SOS_log is not None:
        if SOS_log[ext].data is not None:
            try:
                val = SOS_log[ext].data[np.where(SOS_log[ext].data['FILENAME'] == ptt.basename(ffile))[0]][col]
            except:
                val = []
        else:
            val = []
    else:
        val = []
    if len(val) == 1:
        return(val[0])
    elif len(val) == 0:
        return('-')
    else:
        return(val[-1])



quals = {'excellent':  1,
         'marginal':  .5,
         'bad':        0,
         'test':      -1}

def update_hdr(mjd,obs,hdr):
    quality = None
    fix_file = ptt.join(os.getenv('SDHDRFIX_DIR'),obs.lower(),'sdHdrfix','sdHdrFix-'+str(mjd)+'.par')
    if ptt.exists(fix_file):
        fix = yanny.read_table_yanny(fix_file, 'OPHDRFIX')
        file_root = 'sdR-??-'+str(hdr['EXPOSURE']).zfill(8)

        updates = fix[fix['fileroot'] == file_root]
        for row in updates:
            if row['keyword'] == 'quality':
                qaulity = row['value']
            hdr[row['keyword']] = row['value']
        file_root = file_root.replace('??', hdr['CAMERAS'].strip())
        updates = fix[fix['fileroot'] == file_root]
        for row in updates:
            if row['keyword'] == 'quality':
                qaulity = row['value']
            hdr[row['keyword']] = row['value']
    return(hdr, quality)

def log_exp(ffile, arc, temp, ref, SOS_log, sos_dir,mjd, obs, long_log = False, new_ref = False, hdrfix = None):
    try:
        hdr = fits.getheader(ffile)
    except:
        try:
            hdr = fits.getheader(ffile.replace('.gz',''))
            ffile = ffile.replace('.gz','')
        except:
            try:
                hdr = fits.getheader(ffile+'.gz')
                ffile = ffile+'.gz'
            except:
                return(None)

    ccddef = 'b1' if obs.lower() == 'apo' else 'b2'

    p98 = sn2 = sn2v2 = sn2m15 = sig = f_rat = w_shift = quality = sky = '-'

    hdr,qual = update_hdr(mjd, obs, hdr)
    flav = hdr.get('FLAVOR','??').lower()
    if flav == 'arc':
        ext = 3
        col = 'WSIGMA'
        sig = get_SOS_log_val(SOS_log, ffile, ext, col)
        if type(sig) != str:
            sig = "{:.2f}".format(sig)
        wframe = ptt.join(sos_dir, str(mjd),'wset-'+str(mjd)+'-'+str(hdr.get('FIELDID','??????')).zfill(6)+'-'+str(hdr.get('EXPOSURE','????????')).zfill(8)+'-'+hdr.get('CAMERAS',ccddef).strip()+'.fits')
        if ptt.exists(wframe):
            frame = fits.getdata(wframe,2)
            if not new_ref:
                try:
                    f_rat = (np.median(frame)/hdr['EXPTIME'])/(ref.arc.f_ref[hdr.get('CAMERAS',ccddef).strip()]/ref.arc.exptime)
                except:
                    f_rat = 0.0
                f_rat = "{:.1f}".format(f_rat)
            else:
                f_rat = (np.median(frame))
        if not new_ref:
            wmid = get_SOS_log_val(SOS_log, ffile, ext, 'WAVEMID')
            if type(wmid)!= str:
                try:
                    w_shift = (wmid - ref.arc.WAVEMID_ref[hdr.get('CAMERAS',ccddef).strip()])
                except:
                    w_shift = 0
                w_shift = "{:.1f}".format(w_shift)
        else:
            w_shift = (get_SOS_log_val(SOS_log, ffile, ext, 'WAVEMID'))
    elif flav == 'flat':
        ext = 2
        col = 'XSIGMA'
        sig = get_SOS_log_val(SOS_log, ffile, ext, col)
        if type(sig) != str:
            sig = "{:.2f}".format(sig)
        try:
            frame = fits.getdata(ptt.join(sos_dir, str(mjd),
                      'tset-'+str(mjd)+'-'+str(hdr['FIELDID']).zfill(6)+'-'+str(hdr.get('EXPOSURE','????????')).zfill(8)+'-'+hdr.get('CAMERAS',ccddef).strip()+'.fits'),0)
            if not new_ref:
                f_rat = (np.median(frame)/hdr['EXPTIME'])/(ref.flat.f_ref[hdr.get('CAMERAS',ccddef).strip()]/ref.flat.exptime)
                f_rat = "{:.1f}".format(f_rat)
            else:
                f_rat = (np.median(frame))
        except:
            pass
    elif (flav == 'bias')  or (flav == 'dark'):
        ext = 1
        col = 'PERCENTILE'
        per = get_SOS_log_val(SOS_log, ffile, ext, col)
        try:
            if type(per) != str:
                p98 = "{:.0f}".format(per[97])
                if not new_ref:
                    if flav == 'dark':
                        f_rat = (per[97]/hdr['EXPTIME'])/(ref.dark.f_ref[hdr.get('CAMERAS',ccddef).strip()]/ref.dark.exptime)
                    else:
                        f_rat = (per[97])/(ref.bias.f_ref[hdr.get('CAMERAS',ccddef).strip()])
                    f_rat = "{:.1f}".format(f_rat)
                else:
                    f_rat = per[97]
        except:
            f_rat = 0
    elif (flav == 'science'):
        ext = 4
        col = 'SN2'
        sn2 = get_SOS_log_val(SOS_log, ffile, ext, col)
        if type(sn2) != str:
            sn2 = "{:.1f}".format(sn2)
        col = 'SN2_V2'
        sn2v2 = get_SOS_log_val(SOS_log, ffile, ext, col)
        if type(sn2v2) != str:
            sn2v2 = "{:.1f}".format(sn2v2)
        col = 'SN2_15'
        sn2m15 = get_SOS_log_val(SOS_log, ffile, ext, col)
        if type(sn2m15) != str:
            sn2m15 = "{:.1f}".format(sn2m15)
            
        sky = get_SOS_log_val(SOS_log, ffile, ext,'SKYPERSEC')
        if type(sky) != str:
            sky = "{:.2f}".format(sky).lstrip('0')
    else:
        ext = None
    if ext is not None:
        quality = get_SOS_log_val(SOS_log,ffile, ext, 'QUALITY')

    exp = pd.Series({'Time':time.ctime(ptt.getctime(ffile)), 'Filename': ptt.basename(ffile),
                     'EXPTIME':hdr.get('EXPTIME',0), 'CCD':hdr.get('CAMERAS','??').strip(), 'EXPID':hdr.get('EXPOSURE','????????'), 'DESIGNID':hdr.get('DESIGNID','????'),
                     'Field': hdr.get('FIELDID','??????'), 'configID':hdr.get('CONFID','??????'), 'flavor': flav, 'Hrt':hdr.get('HARTMANN','?'),
                     'FF':hdr.get('FF','? ? ? ?'), 'FFS':hdr.get('FFS','? ? ? ? ? ? ? ?'), 'NE': hdr.get('NE','? ? ? ?'), arc: hdr.get(arc,'? ? ? ?'),
                     'colA': hdr.get('COLLA',0), 'colB': hdr.get('COLLB',0),'colC': hdr.get('COLLC',0), '98%':p98, 'fratio':f_rat,
                     'W(X)':sig, 'w_shift':w_shift, 'SN2':sn2,'SN2_V2':sn2v2,'SN2_15':sn2m15, 'sky/s':sky,'QUALITY':quality, 'temp':hdr.get(temp,'?')})

    exp['full_time'] = time.ctime(ptt.getctime(ffile))
    if (exp['Hrt'].lower() != 'out') & (exp['Hrt'] != 'Closed, Closed') & (exp['Hrt'] != '?') & (exp['Hrt'] is not None):
        exp['flavor'] = 'hart'
        exp['fratio'] = '-'
    if exp['Hrt'].lower() == 'right':
        exp['Hrt'] = 'R'
    elif exp['Hrt'].lower() == 'left':
        exp['Hrt'] = 'L'
    if exp['Hrt'] == 'Closed, Closed':
        exp['Hrt'] = 'Closed'
    for key in ['DESIGNID', 'Field','configID']:
        if exp[key] == -999:
            exp[key] = '-'
    for key in [arc, 'NE', 'FF', 'FFS']:
        if exp[key] == 0:
            exp[key] = '-'
    if not long_log:
        exp['Filename'] = exp['Filename'].replace('sdR-','').replace('.fit.gz','')
        for key in [arc, 'NE', 'FF', 'FFS']:
            test = np.asarray(exp[key].split())
            failed = np.where(np.logical_or((np.char.upper(test) == 'X'), (np.char.upper(test) == '?')))[0]
            test[failed] = 2
            test = test.astype(int)
            test[failed] = 10
            exp[key] = np.sum(test)
            if ((exp['flavor'].lower() == 'arc') or (exp['flavor'].lower() == 'hart')) & (key in [arc,'NE']):
                continue
            if (exp['flavor'].lower() == 'hart') & (key == 'FFS'):
                continue
            if (exp['flavor'].lower() == 'flat') & (key in ['FF','FFS']):
                continue
            if exp[key] == 0:
                exp[key] = '-'
        exp['lamps'] = ','.join(exp[['FF','NE',arc]].values.astype(str))
        exp = exp.drop(['FF','NE',arc])
        if exp['lamps'] == '-,-,-':
            exp['lamps'] = '  -  '
        exp = exp.rename({'EXPTIME':'EXP', 'flavor':'Flav', 'DESIGNID':'DESIGN'})
        exp['Time'] = exp['Time'].split()[-2]
        exp['cols(A,B,C)'] = ','.join(exp[['colA','colB','colC']].values.astype(str))
        exp = exp.drop(['colA','colB','colC'])
        if exp['QUALITY'] != '-':
            try:
               exp['QUALITY'] = quals[exp['QUALITY']]
               if exp['QUALITY'] == 0.5:
                   exp['QUALITY'] = '.5'
            except:
               exp['QUALITY'] = '?'
        exp = exp.rename({'QUALITY':'Q'})
        if type(exp['EXP']) != str:
           exp['EXP'] = int(exp['EXP'])
    exp = exp[empty_log(arc, long_log).columns]
    return(exp)


def empty_log(arc, long_log = False):
    if long_log:
        cols = ['full_time','Time','Filename','CCD','EXPID','DESIGNID','Field','configID','EXPTIME','flavor','Hrt','FF','FFS','NE',arc,'colA','colB','colC','temp','98%','fratio','W(X)','w_shift','sky/s','SN2','SN2_V2','SN2_15','QUALITY']
    else:
        cols = ['full_time','Time','Filename','CCD','EXPID','Q','DESIGN','Field','configID','cols(A,B,C)','temp','EXP','Flav','Hrt','FFS','lamps','98%','fratio','W(X)','w_shift','sky/s','SN2','SN2_V2','SN2_15']
    return(pd.DataFrame(columns=cols,index=[-1]))


def merge_ccd(log, col, ccds):
    blue = log[log.CCD == ccds[0]].copy()
    if len(blue) == 0:
        blue = '-'
    else:
        blue = (blue.iloc[0])[col]
    red = log[log.CCD == ccds[1]].copy()
    if len(red) == 0:
        red = '-'
    else:
        red = (red.iloc[0])[col]
    return(','.join([str(blue),str(red)]))


def get_hartmann_logs(hart_time, obs='lco'):
    logdir = '/data/logs/tron/hartmann'
    hlogs = []

    hart_time = datetime.datetime.strptime(hart_time, "%a %b %d %H:%M:%S %Y")
    test_logs = []
    test_logs.append((hart_time+datetime.timedelta(days=1)).strftime("%Y-%m-%d")+'*')
    test_logs.append((hart_time).strftime("%Y-%m-%d")+'*')
    test_logs.append((hart_time+datetime.timedelta(days=-1)).strftime("%Y-%m-%d")+'*')

    if obs == 'lco':
        cams = ['b2','r2']
        spec = 'sp2'
    else:
        cams = ['b1','r1']
        spec = 'sp1'

    for tl in test_logs:
       for tll in glob(ptt.join(logdir,tl)):
            v2 = False
            with open(tll) as tlog:
                lines = tlog.readlines()
                for i, line in enumerate(lines):
                    line_time = line.split(',')[0]
                    if line_time == '\n':
                        continue
                    line_time = datetime.datetime.strptime(line_time, "%Y-%m-%d %H:%M:%S")
                    if abs((line_time - hart_time).total_seconds()) < 600.0:
                        if (((cams[0] in line) or (cams[1] in line) or (spec in line) or ('collimator' in line) or
                            (('>' in line) and ('collimate' in line)) or ('status="moving"' in line)) and
                            (('text' not in line) and (f'specs={spec}' not in line) and (f'cameras={cams[0]},{cams[1]}' not in line))):
                            
                            if ('status="moving"' in line):
                                v2 = True
                                hlogs.append(line)
                            elif v2:
                                hlogs.append(line.replace('\n',''))
                            else:
                                hlogs.append(line.replace('\n',''))
                        if not v2:
                            if ('Not moving collimator' in line):
                                hlogs.append(line.replace('\n',''))
                            if ('Adjusting collimator' in line):
                                hlogs.append(line.replace('\n',''))
                            if ('Hartmann collimation failed' in line):
                                hlogs.append(line)
                            if (('status=idle' in line) and ('Adjusting collimator' in lines[i-1])):
                                hlogs[-1] += '\n'
                            if ('error="Found error' in line):
                                hlogs.append(line.replace('\n',''))

    return(hlogs)


def parse_hartmann_logs(hlogs, exps, long_log = False):
    if long_log:
        cols = {'flav':'flavor', 'exp':'EXPID', 'design': 'DESIGNID', 'configID':'configID', 'Field':'Field'}
    else:
        cols = {'flav':'Flav', 'exp':'Filename', 'design': 'DESIGN', 'configID':'configID','Field':'Field'}
    exps = exps[exps[cols['flav']] == 'hart']
    exps = exps.reset_index()
    idx = []
    for i, exp in exps.iterrows():
        if i == 0:
            continue
        if exp.Hrt == 'R':
            if (int(exp.EXPID) -1 in exps.EXPID.values.astype(int)) & (exps.iloc[i-1].Hrt == 'L'):
                idx.append(exp.name)
    exps = exps.drop(index=idx)

    harttime = []
    for i, exp in exps.iterrows():
        harttime.append(datetime.datetime.strptime(exp.full_time.split()[-2], "%H:%M:%S"))
    harttime = np.asarray(harttime)

    hartmann_logs = pd.DataFrame()
    if len(hlogs) == 0:
        return(hartmann_logs)
    hart = None
    for i, le in enumerate(hlogs):
        var = le.split('=')[0].split()[-1]
        if '>' in le:
            if i > 0 and hart is not None:
               if 'status' in hart.index:
                   hart.status = ','.join(hart.status)
               hartmann_logs = pd.concat([hartmann_logs, pd.DataFrame([hart])], ignore_index = True)

            hart_time = datetime.datetime.strptime(le.split(',')[0].split()[1], "%H:%M:%S")
            delta = np.abs(harttime - hart_time)
            match = exps.iloc[np.argmin(delta)].copy()
            hart = pd.Series({'time': le.split(',')[0].split(' ')[1], 'exp0':match[cols['exp']],'Field':match[cols['Field']], 'DESIGN':match[cols['design']],
                              'configID':match[cols['configID']], 'flag':le.split('collimate')[-1], 'temp['+chr(176)+'C]':match['temp']})
        elif ('MeanOffset' in le):
            ccd = var.split('MeanOffset')[0]
            hart[ccd.upper()+'off'] = (le.split('=')[1]).split(',')[0]
            hart[ccd.upper()] = 'out' if 'Out of focus' in le else 'in'
        elif ('PistonMove' in le):
            ccd = var.split('PistonMove')[0]
            hart[ccd.upper()+'Pmove'] = le.split('=')[1]
        elif ('RingMove' in le):
            ccd = var.split('RingMove')[0]
            hart[ccd.upper()+'Rmove'] = le.split('=')[1]
        elif ('Residuals' in le):
            out = le.split('=')[1]
            out = ','.join([ x for x in out.split(',') if x.lstrip('-').replace('.','',1).isdigit()])
            hart['Residuals'] = out
           
        elif ('AverageMove' in le):
            hart['AvgMove'] = le.split('=')[1]
        elif ('error' in le):
            if 'status' not in hart.index:
                hart['status'] = []
            hart['status'].append(le.split('=')[-1].replace('"', '').replace('\n',''))
            hart['status'] = [*dict.fromkeys(hart['status'])]
        elif ('text' in le):
            if 'status' not in hart.index:
                hart['status'] = []
            stat = le.split('=')[1].replace('"','').replace('\n','')
            if ('ignore_residuals' in stat) or ('Not moving' in stat) :
                continue
            hart['status'].append(stat)
            hart['status'] = [*dict.fromkeys(hart['status'])]
        else:
            if hart is None:
                continue
            if 'status' not in hart.index:
                hart['status'] = []
            stat = le.split('=')[1].replace('"','').replace('\n','')
            hart['status'].append(stat)
            hart['status'] = [*dict.fromkeys(hart['status'])]
    if 'status' in hart.index:
        hart.status = ','.join(hart.status)
    hartmann_logs = pd.concat([hartmann_logs, pd.DataFrame([hart])], ignore_index = True)
    hartmann_logs = hartmann_logs.fillna('')
    cols_at_end = ['status']
    hartmann_logs = hartmann_logs[[c for c in hartmann_logs if c not in cols_at_end]
                                + [c for c in cols_at_end if c in hartmann_logs]]
    return(hartmann_logs)


def print_hart(log, obs, hart_table, long_log=False):
    hart_logs = []
    Collimation = pd.DataFrame()
    fl = 'flavor' if (long_log is True) else 'Flav'
    for i, exp in log[log[fl] == 'hart'].iterrows():
        hart_logs.extend(get_hartmann_logs(exp.full_time, obs=obs))
    hart_logs = [*dict.fromkeys(hart_logs)]
    if len(hart_logs) == 0:
        hart_table = False
    if hart_table:
        hart_logs_tab = parse_hartmann_logs(hart_logs, log, long_log)
    print('    ---- BOSS Spectrographs Collimation ----')
    if hart_table:
        hart_logs_tab = hart_logs_tab.replace({'\"':''}, regex=True)
        head = hart_logs_tab.to_string(index=False).split('\n')
        if len(head[0])> 200:
            hart_table=False
        else:
            if len(hart_logs_tab) == 0:
                hlog = '\n'
            else:
                hlog = '\n    '.join(head[1:])
                hlog = '    '+hlog
            head = head[0]
            print('    '+'- '*int(len(head)/2+1))
            print('    '+head)
            print('    '+'- '*int(len(head)/2+1))
            print(hlog)
    if not hart_table:
        for hl in hart_logs:
            print('    '+hl)
    print('\n')

def print_SOSwarn(sos_dir, mjd):
    print('\n    ---- SOS ERRORS/WARNGINGS ----')
    try:
        messages = fits.getdata(ptt.join(sos_dir,str(mjd),'logfile-'+str(mjd)+'.fits'), 5)['TEXT']
        for m in messages:
            if 'WARNING' in m:
                ms = m.split('WARNING:')
                print('      '+ms[0]+colored( 'WARNING:','yellow')+ms[1])
            elif 'ABORT' in m:
                ms = m.split('ABORT:')
                print('      '+ms[0]+colored( 'ABORT:', 'red')+ms[1])
            else:
                print('      '+m)
        print('\n')
    except:
        print('\n')


def print_summary(Datadir, sos_dir, mjd, vers2d, run2d, log, quals, obs, arc, ref, long_log=False):
    log = log.drop(columns = ['EXPID', 'full_time'])
    log = log.replace('science','sci')
    log = log.rename(columns = {'temp':'T['+chr(176)+'C]'})

    print('---- BOSS Data Summary ----')
    print('Reading FITS from: '+ptt.join(Datadir, mjd)+',  '+ptt.join(sos_dir, mjd))
    print('Using '+ptt.join(sos_dir,str(mjd),'logfile-'+str(mjd)+'.fits')+ ' from IDLSPEC2D version ' +vers2d +'  (RUN2D '+run2d+')')
    head = log.to_string().split('\n')
    if log.index[0] == -1:
        log = '\n'
    else:
        log  = '\n'.join(head[1:])
    head =head[0]
    print('- '*int(len(head)/2+1))
    print('n'+head[1:])
    print('- '*int(len(head)/2+1))
    print(log)

    if (log != '\n'):
        print('- '*int(len(head)/2+1))
        ntt = 'MC1TEMDN: sp1 mech Median temp' if obs == 'apo' else "COLLT: Collimator temperature"
        if (not long_log):
            print("      Q = Quality "+str(quals)+";    T["+chr(176)+"C] = " +ntt)
            nt = ' (pseudo at LCO)' if obs == 'lco' else ''
            print("      FFS = Flat Screens"+nt+" {'open': 0, 'closed': 1, 'Failure - X(?)': 10};    lamps = 4,4,4 (FF, NE, "+arc+") {'on': 1,'off': 0,'Failure - X(?)': 10}")
        else:
            print("      T["+chr(176)+"C] = " +ntt)
        print("      fratio & w_shift(Anstrom of central pixel) are comparison to reference data: flat-"+str(ref.flat.mjd)+" arc-"+str(ref.arc.mjd)+" bias-"+str(ref.bias.mjd)+" dark-"+str(ref.dark.mjd))


def built_short_log(log, ccds ):
    log_short = pd.DataFrame()
    for eid in log.EXPID.unique():
        tlog = log[log.EXPID == eid].iloc[0].copy()
        tlog['Filename'] = tlog['Filename'].split('-')[1].split('.')[0]
        flav = tlog['Flav']
        if (flav == 'arc') or (flav == 'flat'):
             tlog['W(X)'] = merge_ccd(log[log.EXPID == eid], 'W(X)', ccds)
             tlog['fratio'] = merge_ccd(log[log.EXPID == eid], 'fratio', ccds)
             if flav == 'arc':
                 tlog['w_shift'] = merge_ccd(log[log.EXPID == eid], 'w_shift', ccds)
        elif flav == 'science':
             tlog['SN2'] = merge_ccd(log[log.EXPID == eid], 'SN2', ccds)
             tlog['SN2_V2'] = merge_ccd(log[log.EXPID == eid], 'SN2_V2', ccds)
             tlog['SN2_15'] = merge_ccd(log[log.EXPID == eid], 'SN2_15', ccds)
             tlog['sky/s'] = merge_ccd(log[log.EXPID == eid], 'sky/s', ccds)
        elif (flav == 'bias') or (flav == 'dark'):
             tlog['98%'] = merge_ccd(log[log.EXPID == eid], '98%', ccds)
             tlog['fratio'] = merge_ccd(log[log.EXPID == eid], 'fratio', ccds)
        for col in ['CCD']:
             tlog[col] = merge_ccd(log[log.EXPID == eid], col, ccds)
        log_short = pd.concat([log_short,pd.DataFrame([tlog])], ignore_index=True)
    log = log_short
    log = log.rename(columns={'W(X)':'W(X)SIG'})
    return(log)

def get_run2d(sos_log):
    try:
        vers2d = sos_log[0].header['VERS2D']
    except:
        vers2d = ''
    try:
        run2d = sos_log[0].header['RUN2D']
    except:
        run2d = vers2d
    if run2d == '':
        run2d = vers2d
    if run2d == '':
        run2d = os.getenv('IDLSPEC2D_VER',default = '')
        if run2d != '':
            vers2d = os.popen('cd $IDLSPEC2D_DIR ; git describe --tags --abbrev=0').read().rstrip("\n")

    return(vers2d, run2d)


def send_email(obs, mjd, raw_output, email):
    try:
        conv = Ansi2HTMLConverter(dark_bg=False, scheme='xterm')
        html_body = conv.convert(raw_output, full=True)
        ANSI_ESCAPE_RE = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
        
        # Original HTML from ansi2html
        html_output = conv.convert(raw_output, full=True)
        soup = BeautifulSoup(html_output, "html.parser")

        # Remove background color rules
        for cls in soup.select("[class*='background']"):
            cls['class'] = [c for c in cls.get('class', []) if 'background' not in c]

        # Remove background-color CSS rules
        if soup.style:
            style_lines = soup.style.string.splitlines()
            filtered = [line for line in style_lines if "background-color" not in line]
            soup.style.string = "\n".join(filtered)
        html_body = str(soup)
        
    except:
        ANSI_ESCAPE_RE = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')

        html_body = ANSI_ESCAPE_RE.sub('', raw_output)
        
        
    if obs.upper() == "LCO":
        sender = "sdss-alerts@lco.cl"
        client = "smtp-02.lco.cl:25"
    else:
        sender = "sdss5-bhm@apo.nmsu.edu"
        client = "localhost" #"mail.apo.nmsu.edu"
    
    msg = MIMEMultipart("alternative")
    msg["Subject"] = f"{obs.upper()} SOS BOSS Log MJD:{mjd}"
    msg["From"] = sender
    msg["To"] = email

    # Add HTML and/or plain version
    msg.attach(MIMEText(ANSI_ESCAPE_RE.sub('', raw_output), "plain"))
    msg.attach(MIMEText(html_body, "html"))

    # Send email
    with smtplib.SMTP(client) as server:
        server.send_message(msg)



def build_log(mjd, obs, Datadir='/data/spectro/', sos_dir = '/data/boss/sos/', long_log = False, new_ref = False,
              hart=False, hart_table=False, hide_error=False, hide_summary=False, email=None):
    log = pd.DataFrame()
    run2d = vers2d = ''

    if obs == 'apo':
        key  = '*-b1*'
        key2 = '*-r1*'
        arc = 'HGCD'
        temp = 'MC1TEMDN'
        ccds = ['b1','r1']
        ref = apo_ref
    else:
        key  = '*-b2*'
        key2 = '*-r2*'
        arc = 'HEAR'
        temp = 'COLLT'
        ccds = ['b2','r2']
        ref = lco_ref

    rkey = key.replace('*','')
    key2 = key2.replace('*','')

    if ptt.exists(ptt.join(sos_dir,str(mjd),'logfile-'+str(mjd)+'.fits')):
        sos_log = fits.open(ptt.join(sos_dir,str(mjd),'logfile-'+str(mjd)+'.fits'))
    else:
        sos_log = None

    ffiles   = sorted(glob(ptt.join(Datadir,mjd,key)), key=get_key)
    ffiles_num = [get_key(x) for x in ffiles]
    ffiles_2 = sorted(glob(ptt.join(Datadir,mjd,key2)),key=get_key)
    ffiles_num_2 = [get_key(x) for x in ffiles_2]
    missing = np.setdiff1d(ffiles_num_2, ffiles_num)
    for mi in missing:
        ffiles.append(ptt.join(Datadir,mjd,'sdR'+rkey+str(mi).zfill(8)+'.fit'))

    for ffile in sorted(ffiles, key=get_key):
        exp = log_exp(ffile,arc,temp,ref, sos_log,sos_dir,mjd,obs, long_log = long_log, new_ref = new_ref)
        if exp is not None:
            log = pd.concat([log,pd.DataFrame([exp])], ignore_index=True)
        exp = log_exp(ffile.replace(rkey,key2),arc,temp, ref, sos_log,sos_dir,mjd, obs,  long_log = long_log, new_ref = new_ref)
        if exp is not None:
            log = pd.concat([log,pd.DataFrame([exp])], ignore_index=True)

    vers2d, run2d =  get_run2d(sos_log)

    if sos_log is not None:
        sos_log.close()
    if not long_log and len(log) >0:
        log = built_short_log(log, ccds )
    if len(log) == 0:
        log = empty_log(arc, long_log = long_log)
        
        
    if email:
        output_lines = []
        def capture_print(*args, **kwargs):
            s = " ".join(str(arg) for arg in args)
            output_lines.append(s)
        orig_print = builtins.print
        builtins.print = capture_print
        
    try:
        if not hide_summary:
            print_summary(Datadir, sos_dir, mjd, vers2d, run2d, log, quals, obs, arc, ref, long_log=long_log)
        if not hide_error:
            print_SOSwarn(sos_dir, mjd)
        if hart:
            print_hart(log,obs, hart_table, long_log=long_log)
    finally:
        if email:
            builtins.print = orig_print
            raw_output = "\n".join(output_lines)
        else:
            pass

    if email:
        send_email(obs, mjd, raw_output, email)
