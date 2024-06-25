#!/usr/bin/env python3

import pandas as pd
import argparse
try:
    from field import field_to_string as f2s
except:
    def f2s(field):
        return(str(field).zfill(6))
import os.path as ptt
import glob
import time
import smtplib
#from email.message import EmailMessage
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from pydl.pydlutils.yanny import yanny, read_table_yanny
import re
from os import getenv, makedirs, rename
import numpy as np
try:
    import json
    read_json = True
except:
    read_json = False
import datetime
import astropy.time
from collections import OrderedDict

    
jdate = int(float(astropy.time.Time(datetime.datetime.utcnow()).jd) - 2400000.5)

    
try:
    daily_dir = getenv('DAILY_DIR')
except:
    daily_dir = ptt.join(getenv('HOME'),'daily')

def chpc2html(fpath):
    return(fpath.replace("/uufs/chpc.utah.edu/common/home/sdss50","https://data.sdss5.org/sas/"))


def get_nextmjd(run2d, obs, nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par')):
    mod = 'bhm/'+run2d
    try:
        nextmjds = yanny(nextmjd_file)
    except:
        nextmjds = {}
        nextmjds["NEXTMJD"] = Table(names=('module', 'mjd', 'obs'), dtype=('S30', int, 'S3'))
    obss  = np.char.upper(nextmjds["NEXTMJD"]['obs'].astype(str))
    mods = np.char.lower(nextmjds["NEXTMJD"]['module'].astype(str))
    indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]
    if len(indx) == 0:
        mod = 'work/'+run2d
        indx = np.where((obss == obs.upper()) & (mods == mod.lower()))[0]

    if len(indx) == 0:
        nextmjd = 0
    else:
        nextmjd = nextmjds["NEXTMJD"]['mjd'][indx][0]
    return(int(nextmjd))

class Crash_log:
    def __init__(self, step, error,msg=None,line=None, color='DarkOrange'):
        self.step  = step
        self.error = re.compile('.*'+error+'.*')
        self.msg = msg
        self.line = line
        self.color = color
    def check(self,i, line, step):
        if self.step is not None:
            if self.step != step:
                return
        if self.line is not None:
            if i > self.line:
                return
        if self.error.match(line):
            if self.msg is None:
                return(line.replace('\n',''))
            else:
                return(self.msg.format(step=step))
   
   
errors = [Crash_log('spDiag2d','LOCATESKYLINES:.*WARNING: Maximum sky-line shift is.*(DISABLING)'),
          Crash_log('spDiag2d','ABORT: Only            0 sky fibers found',
                    msg='No Sky Fibers Found', color='red'),
          Crash_log('spDiag2d','ABORT: No good flats (saturated?)', color='red'),
          Crash_log('spDiag2d','SPCALIB: .*: .* paired with no arc', color='red'),
          Crash_log('spDiag2d','ABORT: Reject science as too bright: 25-th-percentile =',
                    msg='Reject Bright Science', color='red'),
          Crash_log('spDiag2d','SKYSUBTRACT:.*: Discarding .*(fractional) of the sky pixels as bad',
                    msg='Failed Sky Subtraction', line = -1, color='red'),
          Crash_log('spDiag2d','FITSPECTRARESOL: .*: Calculating the spectra resolution',
                    msg='Failed FITSPECTRARESOL', line = 1, color='red'),
          Crash_log('spDiag2d','EXTRACT_BUNDLE_IMAGE: .*: sigmasize:',
                    msg='Failure Extracting Exposure', line = 1, color='red'),
          Crash_log('spDiag2d','FITMEANX: .*:',msg='Failure in Sky Line Identification',
                    line = 1, color='red'),
          Crash_log('spDiagcomb','RM_SPFLUX_V5:.*: USING XYFIT', color='red',
                    msg='SpectroPhoto Calibration Failure', line = 1),
          Crash_log('spDiagcomb','RM_SPCOMBINE_V5: ABORT: No exposures with SCORE > 0',
                    msg='No Good Exposures', color='Maroon'),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: Rejected  .* of  .* std stars',
                    msg='Failure Combining Exposures', line = 1, color='red'),
          Crash_log('spDiag1d','ZCOMPUTE: .*',msg='Failure in COMPUTECHI2 for ZFIND',
                    line = 1, color='red'),
          Crash_log('spDiag1d','ZFIND: .*',msg='Failure in COMPUTECHI2 for ZFIND',
                                  line = 1, color='red'),
          Crash_log('run_spTrace','Killed', msg='Failed run_spTrace', color='red'),
          Crash_log('spAll','fieldmerge: EXITING!!', color='red'),
          Crash_log('spSpec_reformat', 'read_spAll: ERROR: Missing .*',
                    msg='Failed spSpec_reformat: missing spAll field', color='red')]


py_err = [Crash_log(None,'exception:',
                     msg='Failed {step}', color='red'),
                     Crash_log('spAll','fieldmerge: No valid spAll entries', color='red')]

noerr_cal_b = Crash_log('spDiag2d','SPCALIB: b.*: .* paired with arc',
                        msg='No error', color='green')
noerr_cal_r = Crash_log('spDiag2d','SPCALIB: r.*: .* paired with arc',
                        msg='No error', color='green')

def parse_log(file, custom=None):
    if custom is None:
        complete= {'spfibermap':'Successful completion of readfibermaps',
                   'spDiag2d':'Successful completion of SPREDUCE2D',
                   'spDiagcomb':'Successful completion of SPCOMBINE',
                   'spDiag1d':'Successful completion of SPREDUCE1D',
                   'spXCSAO':'CPU time to compute RVs',
                   'fieldlist':'Successful completion of fieldlist',
                   'spAll':'Successful completion of fieldmerge',
                   'spSpec_reformat':'Successful completion of spSpec_reformat',
                   'spCalib_QA':'SpectroPhoto QA Complete',
                   'run_spTrace':'Successful completion of boss_arcs_to_trace'}#,
    else:
        complete = {'spDiagcomb':'SPSPEC_TARGET_MERGE: Successful completion of spspec_target_merge',
                    'spDiag1d':'Successful completion of SPREDUCE1D',
                    'spXCSAO':'CPU time to compute RVs',
                    'spAll':'Successful completion of fieldmerge',
                    'spSpec_reformat':'Successful completion of spSpec_reformat'}
    for key in complete:
        if key not in file:
            continue
        line = -2 if key == 'spfibermap' else -1
        if not ptt.exists(file):
            return('Gold', f'<b>{key}</b> not yet run')
        with open(file) as f:
            try:
                last_line = f.readlines()[line]
            except:
                last_line = ''
        if complete[key] in last_line:
            with open(file) as f:
                lines = f.readlines()
                lines.reverse()
                
                for i, line in enumerate(lines):
                    for err in errors:
                        msg = err.check(i,line,key)
                        if msg is not None:
                            if key == 'spDiag2d':
                                if 'SPCALIB: b' in msg:
                                    for li, ll in enumerate(lines):
                                        noerr = noerr_cal_b.check(li, ll, key)
                                        if noerr == 'No error':
                                            msg = None
                                elif 'SPCALIB: r' in msg:
                                    for li, ll in enumerate(lines):
                                        noerr = noerr_cal_b.check(li, ll, key)
                                        if noerr == 'No error':
                                            msg = None
                            if msg is not None:
                                return(err.color,msg)
            return('green', None)
        elif key == 'run_spTrace' and '.e.log' in file:
            with open(file) as f:
                lines = f.readlines()
                lines.reverse()
                
                for i, line in enumerate(lines):
                    for err in errors:
                        msg = err.check(i,line,key)
                        if msg is not None:
                            return('red',msg)
            return('green',None)
        else:
            if ((key in ['spfibermap','spXCSAO','fieldlist','spAll','run_spTrace'])
                    and '.log' in file):
                with open(file) as f:
                    lines = f.readlines()
                    lines.reverse()
                    
                    for i, line in enumerate(lines):
                        for err in py_err:
                            msg = err.check(i,line,key)
                            if msg is not None:
                                return(err.color,msg)
                return('Gold', None)
            else:
                if time.time() - ptt.getmtime(file) > 300:
                    # check if log has been updated in last 5 minutes
                    # if not then check for errors
                    with open(file) as f:
                        lines = f.readlines()
                        lines.reverse()
                        for i, line in enumerate(lines):
                            for err in errors:
                                msg = err.check(i,line,key)
                                if msg is not None:
                                    return(err.color,msg)
            return('Gold', None)
    return('Gold',None)
class LogCheck:
    def __init__(self, topdir, run2d, run1d, field, mjd, dither='F',
                 epoch=False, custom = None, mjd1d=None):
        self.topdir = topdir
        self.run2d = run2d
        self.run1d = run1d
        self.field = f2s(field)
        self.mjd = mjd
        self.mjd1d = mjd1d
        self.custom = custom
        self.dither = dither
        self.epoch = epoch

    def html(self, fbase=[], exts =None):
        rs = ''
        note = []
        colors = []
        for i, fb in enumerate(fbase):
            if self.epoch:
                file = ptt.join(self.topdir, self.run2d, self.field, 'epoch',
                                fb.format(field=self.field, mjd=self.mjd))
            elif self.custom is not None:
                file = ptt.join(self.topdir, self.run2d, self.field,
                                fb.format(field=self.field, mjd=self.mjd,
                                          custom=self.custom, mjd1d=self.mjd1d))
            else:
                file = ptt.join(self.topdir, self.run2d, self.field,
                                fb.format(field=self.field, mjd=self.mjd))
            if ptt.splitext(file)[-1] == '.log':
                gf = f'log'
                bf = f"<s style='color:Gold;'>log</s> "
            else:
                if ptt.splitext(file)[-1] == '.gz':
                    ext = ptt.splitext(ptt.splitext(file)[0])[-1]
                else:
                    ext = ptt.splitext(file)[-1]
                if exts is not None:
                    ext = exts[i]
                ext = ext.replace('.','')
                gf = f'{ext}'
                bf = f"<s style='color:Gold;'>{ext}</s> "
            color = 'green'
            
            if ptt.splitext(file)[-1] == '.log':
                color, tnote = parse_log(file, custom=self.custom)
                if 'v6_1' in self.run2d:
                    if 'spCalib_QA' in fb:
                        color = 'green'
                        bf = f"<s style='color:green;'>log</s> "
                        tnote = None
                if tnote is not None:
                    note.append(tnote)
                
            if ptt.exists(file):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file)+f" style='color:{color};'>"+gf+"</A> "
                    colors.append(color)
                    continue
                else:
                    with open(file, 'rb') as ff:
                        if len(ff.readlines()) > 100:
                            rs = rs + "<A HREF="+chpc2html(file)+f" style='color:{color};'>"+gf+"</A> "
                            colors.append(color)
                            continue
            elif '*' in file:
                if len(glob.glob(file)) > 0:
                    rs = rs + "<A HREF="+chpc2html(ptt.dirname(file))+f" style='color:{color};'>"+gf+"</A> "
                    colors.append(color)
                    continue
                    
            if ptt.exists(file.replace('.pdf','.ps')):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file.replace('.pdf','.ps'))+f" style='color:{color};'>"+gf+"</A> "
                    colors.append(color)
                    continue
                elif 'spDiagcomb' in fbase[0]:
                    if 'green' in colors[0]:
                        color = 'green'
                        bf = bf.replace('color:Gold','color:green')
            rs = rs + bf
            colors.append(bf)
        if self.dither == 'T':
            rs = rs.replace('color:red','color:#FF0000')
            rs = rs.replace('color:Gold','color:#FFD700')
        if 'color:red' in rs:
            rs  = rs.replace('color:Gold','color:red').replace('color:green','color:red')

        return(rs, note)


def CheckRedux(topdir, run2d, run1d, field, mjd, obs, dither = 'F', epoch=False,
               plan2d=None, custom = None, mjd1d = None):
    lc = LogCheck(topdir, run2d, run1d, field, mjd, dither = dither,
                  epoch=epoch, custom=custom, mjd1d=mjd1d)
    fmjd = pd.Series({}, dtype=object)
    note = OrderedDict()
    fmjd['Field']  = field
    if custom is not None:
        fmjd['MJD']    = mjd1d
    else:
        fmjd['MJD']    = mjd
    fmjd['OBS']    = obs.upper()
    fmjd['Dither'] = dither
    if custom is None:
        fmjd['redux'], _ = lc.html(['redux-{field}-{mjd}.e','redux-{field}-{mjd}.o'])
    else:
        fmjd['redux'], _ = lc.html(['redux_{field}-{mjd}.e','redux_{field}-{mjd}.o',
                                    'redux_{field}-{mjd}_{mjd1d}.e','redux_{field}-{mjd}_{mjd1d}.o'])

    if epoch:
        plan2d = ['../'+x for x in plan2d]
        exts = ['2d']*len(plan2d)
        exts.append('comb')
        plan2d.append('spPlancombepoch-{field}-{mjd}.par')
        fmjd['plans'], _ = lc.html(plan2d, exts=exts)
#        fmjd['plans'], _ = lc.html(['spPlancombepoch-{field}-{mjd}.par'], exts=['comb'])
        fmjd['spfibermap'] = '-'
        note['spfibermap'] = []
        fmjd['spreduce2D'] = '-'
        note['spreduce2D'] = []
    elif custom is not None:
        fmjd['plans'], _ = lc.html(['spPlanCustom-{custom}-{mjd}.par'])
        

    else:
        fmjd['plans'], _ = lc.html(['spPlan2d-{field}-{mjd}.par',
                                    'spPlancomb-{field}-{mjd}.par'],exts=['2d','comb'])
        fmjd['spfibermap'], note['spfibermap'] = lc.html(['spfibermap-{field}-{mjd}.log',
                                                          'spfibermap-{field}-{mjd}.fits'])
        fmjd['spreduce2D'], note['spreduce2D'] = lc.html(['spDiag2d-{field}-{mjd}.log',
                                                          'spDiag2d-{field}-{mjd}.pdf'])
    if custom is None:
        fmjd['specombine'], note['specombine'] = lc.html(['spDiagcomb-{field}-{mjd}.log',
                                                          'spDiagcomb-{field}-{mjd}.pdf',
                                                          'spFluxdistort-{field}-{mjd}.pdf',
                                                          'spSN2d-{field}-{mjd}.pdf'])
        fmjd['spreduce1d'], note['spreduce1d'] = lc.html([run1d+'/spDiag1d-{field}-{mjd}.log'])
        fmjd['spXCSAO'],    note['spXCSAO']    = lc.html([run1d+'/spXCSAO-{field}-{mjd}.log'])
        fmjd['Fieldlist'],  note['Fieldlist']  = lc.html(['fieldlist-{field}-{mjd}.log'])
    else:
        fmjd['specombine'], note['specombine'] = lc.html(['spDiagcomb-{field}-{mjd}.log',
                                                          'spSN2d-{field}-{mjd1d}.pdf'])

        fmjd['spreduce1d'], note['spreduce1d'] = lc.html([run1d+'/spDiag1d-{field}-{mjd1d}.log',
                                                          run1d+'/spDiag1d-{field}-{mjd1d}.pdf'])
        fmjd['spXCSAO'],    note['spXCSAO']    = lc.html([run1d+'/spXCSAO-{field}-{mjd1d}.log'])
    
    
    
    if not epoch:
        fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd}.log',
                                            '../spectra/full/{field}/{mjd}/spAll-{field}-{mjd}.fits.gz'])
        fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd}.log',
                                                         f'../../images/{run2d}/{run1d}/'+'{field}-{mjd}',
                                                         f'../spectra/full/{field}/{mjd}'],
                                                         exts=['.log','.png','.fits'])
    if epoch:
        fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd}.log',
                                            '../../epoch/spectra/full/{field}/{mjd}/spAll-{field}-{mjd}.fits.gz'])
        fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd}.log',
                                                         f'../../../images/{run2d}/epoch/{run1d}/'+'{field}-{mjd}',
                                                         f'../../epoch/spectra/full/{field}/{mjd}/*.fits'],
                                                         exts=['.log','.png','.fits'])
    elif custom is not None:
        fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd1d}.log',
                                             '../spectra/full/{field}/{mjd1d}/spAll-{field}-{mjd1d}.fits.gz'])
        fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd1d}.log',
                                                         f'../../images/{run2d}/{run1d}/'+'{field}-{mjd1d}',
                                                         f'../spectra/full/{field}/{mjd1d}'],
                                                         exts=['.log','.png','.fits'])
    else:
        fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd}.log',
                                            '../spectra/full/{field}/{mjd}/spAll-{field}-{mjd}.fits.gz'])
        fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd}.log',
                                                         f'../../images/{run2d}/{run1d}/'+'{field}-{mjd}',
                                                         f'../spectra/full/{field}/{mjd}'],
                                                         exts=['.log','.png','.fits'])
        fmjd['SpCalib'],    note['SpCalib'] = lc.html(['spCalib_QA-'+run2d+'-{field}-{mjd}.log',
                                                       'spCalib_QA-'+run2d+'-{field}-{mjd}.pdf'])

    fmjd['Note'] = []
    nep = False
    for key in note:
        if note[key] is not None:
            if 'No Good Exposures' in ', '.join(note[key]):
                nep = True
            if nep is True:
                fmjd[key] = fmjd[key].replace('color:red','color:Maroon').replace('color:Gold','color:Maroon')
            fmjd['Note'].append(', '.join(note[key]))
    fmjd['Note'] = list(dict.fromkeys(fmjd['Note'])) #list(set(fmjd['Note']))
    try:
        fmjd['Note'].remove('')
    except:
        pass
#    if 'No Good Exposures' in fmjd['Note']:
#        fmjd['Note'] = ['No Good Exposures']
    fmjd['Note'] = ', '.join(fmjd['Note'])
    
    return(fmjd)


def daily_log_html(obs, mjd, topdir=None, run2d=None, run1d=None, redux=None,
                    email=False, epoch = False, custom=None):
    obs = np.atleast_1d(obs).tolist()
    if topdir is None:
        topdir = getenv('BOSS_SPECTRO_REDUX')
    if run2d is None:
        run2d = getenv('RUN2D')
    if run1d is None:
        run1d = getenv('RUN1D')
    if redux is None:
        if epoch:
            redux = glob.glob(ptt.join(topdir, run2d, '??????','epoch', f'spPlancombepoch-??????-{mjd}.par'))
        elif custom is not None:
            if obs[0].lower() == 'lco':
                custom = custom+'_lco'
            redux = glob.glob(ptt.join(topdir, run2d, custom, f'spPlanCustom-{custom}-{mjd}.par'))
        else:
            redux = glob.glob(ptt.join(topdir, run2d, '??????', f'spPlan2d-??????-{mjd}.par'))

    if epoch:
        plans = glob.glob(ptt.join(topdir, run2d,'??????','epoch', f'spPlancombepoch-??????-{mjd}.par'))
    elif custom is not None:
        plans = glob.glob(ptt.join(topdir, run2d, custom, f'spPlanCustom-{custom}-{mjd}.par'))
    else:
        plans = glob.glob(ptt.join(topdir, run2d, '??????', f'spPlan2d-??????-{mjd}.par'))


    html = pd.DataFrame()
    for r in redux:
        if epoch:
            plan = r.replace('redux-','spPlancombepoch-')
        elif custom is not None:
            plan = r.replace('redux-', 'spPlanCustom-')
        else:
            plan = r.replace('redux-','spPlan2d-')#+'.par'
        if '.par' not in plan:
            plan = plan+'.par'
        yplan = yanny(plan)
        hdr = yplan.new_dict_from_pairs()
        if 'OBS' in hdr.keys():
            thisobs = hdr['OBS'].lower()
        else:
            thisobs = 'apo'
        if thisobs not in obs:
            continue
        if epoch:
            plan2d = hdr['planfile2d'].replace("'",'').split()
        else:
            plan2d = None
        if 'DITHER' in hdr.keys():
            dither = hdr['DITHER']
        elif epoch:
            dither = 'F'
        elif custom is not None:
            dither = 'F'
        else:
            dither = '?'
        rs = ptt.basename(r).split('-')
        field = rs[1]
        mj = rs[2].split('.')[0]
        if int(mj) != int(mjd):
            continue
        if custom is not None:
            mjd1ds = np.unique(read_table_yanny(plan,'COADDPLAN')['EPOCH_COMBINE'].data).astype(str)
        else:
            mjd1ds = [None]
        for mjd1d in mjd1ds:
            fhtml = CheckRedux(topdir, run2d, run1d, field, mj, thisobs,
                               dither=dither, epoch=epoch, plan2d=plan2d,
                               custom = custom, mjd1d=mjd1d)
        
            html = pd.concat([html, pd.DataFrame([fhtml])])

    if email:
        body = [f"<h3>RUN2D: {run2d}",
                f"Observatory: {','.join(obs)}", f"MJD: {mjd}"]
    else:
        body = [f"<h3>RUN2D: {run2d}",
                f"Observatory: {','.join(obs)}", f"MJD: {mjd}"]
    for ob in obs:
        if custom is None:
            SOS_log = ptt.abspath(ptt.join(topdir,'..','sos',ob.lower(),f"{mjd}",f"logfile-{mjd}.html"))
            if ptt.exists(SOS_log):
                SOS_log = f"<a HREF={chpc2html(SOS_log)}>Log</a>" if ptt.exists(SOS_log) else "N/A"
            else:
                sos_dir = 'BOSS_SOS_N' if ob.lower() == 'apo' else 'BOSS_SOS_S'
                SOS_log = ptt.abspath(ptt.join(getenv(sos_dir),f"{mjd}",f"logfile-{mjd}.html"))
                SOS_log = f"<a HREF={chpc2html(SOS_log)}>Log</a>" if ptt.exists(SOS_log) else "N/A"
            body.append(f"{ob.upper()} SOS: {SOS_log}")

            transferlog_json = "/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/staging/{obs}/atlogs/{mjd}/{mjd}_status.json"
            nightlogs = "/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/staging/{obs}/reports/mos/{th}"
            if read_json:
                if ptt.exists(transferlog_json.format(obs=ob.lower(), mjd=mjd)):
                    with open(transferlog_json.format(obs=ob.lower(), mjd=mjd)) as lf:
                        nightlog = json.load(lf)['logfile']
                        if nightlog is None:
                            nightlogh =  "<i>Missing</i>"
                        else:
                            nightlog = nightlogs.format(obs=ob.lower(),th = nightlog)
                            nightlogh = f"<a HREF={chpc2html(nightlog)}>{ptt.basename(nightlog)}</a>"
                    body.append(f"{ob.upper()} Night Log: {nightlogh}")
                elif len(plans) > 0:
                    body.append(f"{ob.upper()} Night Log: <i>Missing</i>")
                else:
                    body.append(f"{ob.upper()} Night Log: ???")

            spTrace = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.o.log"))
            color,_ = parse_log(spTrace)
            spTracee = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.e.log"))
            color2,_ = parse_log(spTracee)
            if color2 != 'green':
                color = 'red'
            if run2d == 'master':
                if mjd < 60418:
                    color = 'green'
            if epoch:
                color = 'green'
            spTrace = f"<a HREF={chpc2html(spTrace)} style='color: {color};'>o.log</a>" if ptt.exists(spTrace) else "N/A"

            spTrace1 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.e.log"))
            spTrace1 = f"<a HREF={chpc2html(spTrace1)} style='color: {color};'>e.log</a>" if ptt.exists(spTrace1) else "N/A"

            spTrace2 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"arcs_{mjd}_{ob.lower()}.html"))
            spTrace2 = f"<a HREF={chpc2html(spTrace2)} style='color: {color};'>Plots</a>" if ptt.exists(spTrace2) else "N/A"
            body.append(f"{ob.upper()} spTrace: {spTrace} {spTrace1} {spTrace2}")
    
        else:
            reduxb = ptt.abspath(ptt.join(topdir, run2d, custom, f'redux_{custom}-{mjd}'))
            color, _ = parse_log(reduxb.replace('redux_','spDiagcomb-')+'.log',custom=custom)
            reduxo = reduxb+'.o'
            reduxo = f"<a HREF={chpc2html(reduxo)} style='color: {color};'>o</a>"
            reduxe = reduxb+'.e'
            reduxe = f"<a HREF={chpc2html(reduxe)} style='color: {color};'>e</a>"
            body.append(f"Redux Coadd: {reduxo} {reduxe}")
    # spAll
    if epoch:
        spAll = ptt.join(topdir,run2d,'epoch',f'spAll-{run2d}.fits.gz')
    elif custom is not None:
        spAll = ptt.join(topdir,run2d,f'spAll-{run2d}-{custom}.fits.gz')
    else:
        spAll = ptt.join(topdir,run2d,f'spAll-{run2d}.fits.gz')
    if ptt.exists(spAll):
        if email:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> ({time.ctime(ptt.getmtime(spAll))})"
        else:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> <span id='spall'></span>"
        body.append(spallh)
    if epoch:
        spAll = ptt.join(topdir,run2d,'epoch',f'spAll-lite-{run2d}.fits.gz')
    elif custom is not None:
        spAll = ptt.join(topdir,run2d,f'spAll-lite-{run2d}-{custom}.fits.gz')
    else:
        spAll = ptt.join(topdir,run2d,f'spAll-lite-{run2d}.fits.gz')
    if ptt.exists(spAll):
        if email:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> ({time.ctime(ptt.getmtime(spAll))})"
        else:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> <span id='spall-lite'></span>"
        body.append(spallh)

    # fieldlist
    if custom is None:
        if epoch:
            flist = ptt.join(topdir,run2d,'epoch',f'fieldlist-{run2d}.fits')
        else:
            flist = ptt.join(topdir,run2d,f'fieldlist-{run2d}.fits')
        if ptt.exists(flist):
            if email:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (fits)</a> ({time.ctime(ptt.getmtime(flist))})"
            else:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (fits)</a> <span id='fieldlistfits'></span>"
            body.append(flisth)
 
        if epoch:
            flist = ptt.join(topdir,run2d,'epoch',f'fieldlist.html')
        else:
            flist = ptt.join(topdir,run2d,f'fieldlist.html')
        if ptt.exists(flist):
            if email:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (html)</a> ({time.ctime(ptt.getmtime(flist))})"
            else:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (html)</a> <span id='fieldlisthtml'></span>"
            body.append(flisth)
        
    body[-1] = body[-1]+"</h3>"
    body.append(html.to_html(index=False, escape=False, justify="center").replace('<td>', '<td align="center">'))

    return(body)

def daily_log_to_file(obs, mjd, topdir=None, run2d=None, run1d=None, redux=None,
                      html_log=None, summary=True, epoch=False, custom = None):
    obs = np.atleast_1d(obs).tolist()
    if run2d is None:
        run2d = getenv('RUN2D')
    outdir = ptt.join(daily_dir, 'logs', 'Status', 'daily', run2d)
    if epoch:
        outdir = ptt.join(daily_dir, 'logs', 'Status', 'epoch', run2d)
    elif custom is not None:
        outdir = ptt.join(daily_dir, 'logs', 'Status', custom, run2d)
    makedirs(outdir, exist_ok=True)

    for obs in obs:
        if html_log is None:
            body = daily_log_html(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d,
                                   redux=redux, epoch=epoch, custom = custom)
        else:
            body = html_log
    
        with open(ptt.join(outdir,f'{mjd}-{obs.upper()}.html'), 'w') as f:
            f.write("<script src='https://code.jquery.com/jquery-3.6.0.min.js'></script>\n")
            f.write("<script src='file_status.js'></script>\n")
            f.write("<br>\n".join(body))
            f.write(' <footer>\n')
            f.write('   <hr>\n')
            f.write('   last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                    str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo)+'\n')
            f.write(' </footer>\n')
    
    daily_log_js(outdir, topdir, run2d, epoch = epoch, custom=custom)
    if summary:
        daily_log_index(outdir, run2d, epoch = epoch, custom=custom)
    
def daily_log_js(directory, topdir, run2d, epoch=False, custom=None):
    if topdir is None:
        topdir = getenv('BOSS_SPECTRO_REDUX')

    cmd = []

    cmd.append("function getlastmod(url,id) {")
    cmd.append("    var req = new XMLHttpRequest();")
    cmd.append("    req.open('HEAD',url, true);")
    cmd.append("    req.onreadystatechange = function() {")
    cmd.append("        if(req.readyState === req.HEADERS_RECEIVED) {")
    cmd.append("            console.log(new Date(req.getResponseHeader('Last-Modified')).toString());")
    cmd.append("            document.getElementById(id).innerHTML = new Date(req.getResponseHeader('Last-Modified')).toString();")
    cmd.append("        }")
    cmd.append("    }")
    cmd.append("    req.send();")
    cmd.append("}")
    cmd.append("")
    cmd.append("$(document).ready(function() {")

    for filep in [f'spAll-{run2d}.fits.gz',f'spAll-lite-{run2d}.fits.gz',
                         f'fieldlist-{run2d}.fits',f'fieldlist.html']:
        if epoch:
            filep = ptt.join(topdir, run2d, 'epoch', filep)
        else:
            filep = ptt.join(topdir, run2d, filep)
        if 'spAll-lite' in filep:
            filet = 'spall-lite'
        elif 'spAll' in filep:
            filet = 'spall'
        elif 'fieldlist.html' in filep:
            filet = 'fieldlisthtml'
        else:
            filet = 'fieldlistfits'
            
            
        
        cmd.append(f"    getlastmod('{chpc2html(filep)}', '{filet}');")
    cmd.append("});")

    with open(ptt.join(directory, 'file_status.js'), 'w') as f:
        for cl in cmd:
            f.write(cl+'\n')
    

def daily_log_index(directory, RUN2D, epoch = False, custom=None):
    with open(ptt.join(directory,'index.html.tmp'), 'w') as f:
        f.write('<html>\n')
        f.write(' <head>\n')
        if epoch:
            f.write(f'   <title>Epoch BOSS Pipeline Status: {RUN2D}</title>\n')
        elif custom is not None:
            f.write(f'   <title>{custom} BOSS Pipeline Status: {RUN2D}</title>\n')
        else:
            f.write(f'   <title>Daily BOSS Pipeline Status: {RUN2D}</title>\n')
        f.write('   <style>BODY {font: normal 12pt Helvetica,Arial}</style>\n')
        f.write('   <link rel="shortcut icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" type="image/x-icon">\n')
        f.write(' </head>\n')
        f.write(' <body>\n')
        f.write('   <img src="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" style="height: 78; width: auto;">\n')
        if epoch:
            f.write(f'   <h2>Epoch BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')
        elif custom is not None:
            f.write(f'   <h2>{custom} BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')
        else:
            f.write(f'   <h2>Daily BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')

        f.write(' <hr>\n')
        
        f.write(" <div class='row' style='display: flex; margin-left:-5px;margin-right:-5px;'>")
        f.write("  <div class='column' style='flex: 50%; padding: 5px;'>")
        f.write('    <TABLE BORDER=2; padding=15px 10px id="mytable";>\n')
        logs = glob.glob(ptt.join(directory,'?????-???.html'))
        logs = [ptt.basename(x).split('-')[0] for x in logs]
        logs = np.unique(np.asarray(logs)).tolist()

        nextmjd ={}
        if custom is None:
            nextmjd['APO'] = get_nextmjd(RUN2D, 'APO', nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par'))
            nextmjd['LCO'] = get_nextmjd(RUN2D, 'LCO', nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par'))
        else:
            nextmjd['APO'] = jdate
            nextmjd['LCO'] = jdate

        for mjd in sorted(logs,reverse=True):
            obs = sorted(glob.glob(ptt.join(directory,f'{mjd}-???.html')))
            obs = [ptt.basename(x).split('-')[1].split('.')[0] for x in obs]
            
            obs_str = []
            for ob in ['APO','LCO']:
                if ob in obs:
                    color='green'
                    sos=True
                    sptrace=True
                    redux = False

                    with open(ptt.join(directory,f'{mjd}-{ob}.html')) as fl:
                        for line in fl.readlines():
                            if 'color:red;' in line:
                                color='red'
                                break
                            if 'color:DarkOrange;' in line:
                                color='DarkOrange'
                            if 'color:Gold;' in line:
                                if color not in ['DarkOrange']:
                                    color='Gold'
                            if 'color:Maroon;' in line:
                                if color not in ['DarkOrange','Gold']:
                                    color='Maroon'
                            if '???' in line:
                                color='magenta'
                            if 'SOS: N/A' in line:
                                sos = False
                            if f'{ob} spTrace: N/A N/A N/A' in line:
                                if 'v6_1' not in RUN2D:
                                    sptrace = False
                            if 'spfibermap' in line:
                                redux = True
                    
                    if not redux:
                        if nextmjd[ob] <= int(mjd):
                            color='magenta'
                    transferflag = f"/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/staging/{ob.lower()}/atlogs/{mjd}/transfer-{mjd}.done"
                    if ptt.exists(transferflag):
                        if sos is False and sptrace is False and color !='magenta':
                            color = 'blue'
                    elif int(mjd) >= 60150:
                        color='magenta'
                    obs_str.append(f"<A HREF='{mjd}-{ob}.html' style='color: {color};'>{ob}</A>")
                else:
                    obs_str.append(f"<A HREF='{mjd}-{ob}.html' style='color:red;'><S style='color:red;'>{ob}</S></A>")
            f.write(f'      <TR><TD>{mjd}<TD>{obs_str[0]}<TD>{obs_str[1]}</TR>\n')
        f.write('    </TABLE><br>\n')
        f.write('  </div>\n') # end column
        f.write("  <div class='column' style='flex: 50%; padding: 5px;'>")
        f.write('    <TABLE BORDER=2; padding=15px 10px id="mytable";>\n')
        f.write('      <caption style="white-space: nowrap">Key</caption>\n')
        f.write("      <TR><TD><b>Color</b><TD><b>Meaning</b></TR>\n")
        f.write("      <TR><TD><p style='color: magenta;'>MAGENTA</p><TD>Reduction has yet to be started due to incomplete transfer</TR>\n")
        f.write("      <TR><TD><p style='color: red;'>RED</p><TD>ERROR Stopped Reduction</TR>\n")
        f.write("      <TR><TD><p style='color: Maroon;'>MAROON</p><TD>ERROR Stopped Reduction for NO GOOD EXPOSURES</TR>\n")
        f.write("      <TR><TD><p style='color: DarkOrange;'>ORANGE</p><TD>Pipeline ran with errors/warnings</TR>\n")
        f.write("      <TR><TD><p style='color: Gold;'>YELLOW</p><TD>Pipeline is still running</TR>\n")
        f.write("      <TR><TD><p style='color: blue;'>BLUE</p><TD>No Reductions</TR>\n")
        f.write("      <TR><TD><p style='color: green;'>GREEN</p><TD>No issues</TR>\n")
        f.write('    </TABLE><br>\n')
        f.write('  </div>\n') # end column
        f.write(' </div>\n') # end row

        f.write(' </body>\n')
        f.write(' <footer>\n')
        f.write('   <hr>\n')
        f.write('   last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo)+'\n')
        f.write(' </footer>\n')
        f.write('</html>\n')
        
    try:
        rename(ptt.join(directory,'index.html.tmp'), ptt.join(directory,'index.html'))
    except:
        pass
    
    


def daily_log_email(subject, attachment, logger, obs, mjd,
                    email_file = None, topdir=None, epoch=False,
                    run2d=None, run1d=None, content=None, custom=None,
                    from_domain="chpc.utah.edu",  redux = None):
  
    
    body = daily_log_html(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d,
                         redux=redux, email=True, epoch=epoch, custom=custom)
    
    daily_log_to_file(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d,
                     redux=redux, epoch=epoch, custom=custom)#, html_log=body)
    
    try:
        emails = open(email_file).read().splitlines()
    except:
        emails = []
        logger.info(email_file+' does not exist')

    msg = MIMEMultipart("alternative")
    msg['Subject'] = subject
    msg['From'] = f"BOSS Pipeline <{getenv('USER')}@{from_domain}>"
    msg['BCC'] = ', '.join(emails)



    part1 = MIMEText("An html compatible email view is required to view the full email","plain")
    part2 = MIMEText("<br>\n".join(body), "html")
    msg.attach(part1)
    msg.attach(part2)

    if attachment is not None:
        attachment = np.atleast_1d(attachment)
        
        
        for f in attachment:
            with open(f, "rb") as fil:
                part = MIMEApplication(
                    fil.read(),
                    Name=ptt.basename(f)
                )
            # After the file is closed
            part['Content-Disposition'] = 'attachment; filename="%s"' % ptt.basename(f)
            msg.attach(part)
        
    s = smtplib.SMTP('localhost')
    s.send_message(msg)
    s.quit()
    return(None)


if __name__ == '__main__' :
    """
    Build Daily Status emails/htmls
    """
    parser = argparse.ArgumentParser( description='Build/load BOSS Pipeline Status Pages')


    parser.add_argument('--obs',      type=str, help='Observatory for status update', nargs='+', default=['apo','lco'])
    parser.add_argument('--mjd',      type=int, help = 'Update these MJDs', nargs='*')
    parser.add_argument('--mjdstart', type=int, help = 'Starting MJD')
    parser.add_argument('--mjdend',   type=int, help = 'Ending MJD')
    parser.add_argument('--epoch',    action='store_true', help = 'Run for epoch Coadds')
    parser.add_argument('--custom',   type=str, help='Name of custom Coadd', default=None)

    parser.add_argument('--topdir',   type=str, default = getenv('BOSS_SPECTRO_REDUX'),
            help='Optional override value for the environment variable $BOSS_SPECTRO_REDUX')
    parser.add_argument('--run1d',    type=str, default=getenv('RUN1D'),
            help='Optional override value for the enviro variable $RUN1D')
    parser.add_argument('--run2d',    type=str, default=getenv('RUN2D'),
            help='Optional override value for the enviro variable $RUN2D')

    parser.add_argument('--email', action='store_true', help='Send each mjd status as email')
    parser.add_argument('--fast', action='store_true', help='Skip updating index until end')
    args = parser.parse_args()

    if args.run2d is not None:
        if args.run1d is None:
            args.run1d = args.run2d
    
    if args.mjd is None:
        if args.custom is None:
            if args.mjdend is None and args.mjdstart is None:
                args.mjd = [jdate-1, jdate]
            elif args.mjdend is None and args.mjdstart is not None:
                args.mjd = range(args.mjdstart, jdate+1)
            elif args.mjdstart is None and args.mjdend is not None:
                args.mjd = [args.mjdend]
            else:
                args.mjd = range(args.mjdstart, args.mjdend+1)
        else:
            redux = glob.glob(ptt.join(args.topdir, args.run2d, args.custom, f'redux_{args.custom}-?????'))
            redux.extend(glob.glob(ptt.join(args.topdir, args.run2d, args.custom, f'redux_{args.custom}_lco-?????')))
            args.mjd = [int(x.split('-')[-1]) for x in redux]

    for mjd in args.mjd:
        for obs in args.obs:
            print(mjd, obs)
            if args.email:
                if epoch:
                    subject = f'Epoch Status: {mjd} {obs}'
                else:
                    subject = f'Status: {mjd} {obs}'
                daily_log_email(subject, None, None, obs, mjd, content=None,
                            email_file = ptt.join(getenv('HOME'), 'daily', 'etc','emails'),
                            topdir=args.topdir, run2d=args.run2d, run1d=args.run1d,
                            from_domain="chpc.utah.edu",  redux = None, epoch=args.epoch,
                            custom=args.custom)
            else:
                daily_log_to_file(obs, mjd, topdir=args.topdir, run2d=args.run2d,
                            run1d=args.run1d, redux=None, html_log=None,
                            summary=(not args.fast), epoch=args.epoch,
                            custom = args.custom)
    if args.fast:
        if args.epoch:
            dir_ = 'epoch'
        elif args.custom is not None:
            dir_ = custom
        else:
            dir_ = 'daily'
        daily_log_index(ptt.join(daily_dir, 'logs', 'Status', dir_, args.run2d), args.run2d,
                        epoch=args.epoch, custom=args.custom)

