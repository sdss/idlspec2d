#!/usr/bin/env python3

from boss_drp.field import field_to_string as f2s
from boss_drp.field import field_dir, field_spec_dir, field_png_dir
from boss_drp.utils.hash import check_hash
from boss_drp.utils import jdate
from boss_drp import daily_dir

import pandas as pd
import argparse
import os.path as ptt
import glob
import time
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from pydl.pydlutils.yanny import yanny, read_table_yanny
import re
from os import getenv, makedirs, rename, symlink
import numpy as np
try:
    import json
    read_json = True
except:
    read_json = False
import datetime
import astropy.time
from collections import OrderedDict
from bs4 import BeautifulSoup
    

def chpc2html(fpath):
    sas_base_dir = getenv('SAS_BASE_DIR', default=None)
    if sas_base_dir is None:
        sas_base_dir = "/uufs/chpc.utah.edu/common/home/sdss50"
    return(fpath.replace(sas_base_dir,"https://data.sdss5.org/sas/"))


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

class Flag:
    def __init__(self,name,color,code,desc):
        self.name=name
        self.color=color
        self.code=code
        self.desc=desc
    def key(self):
        return(f"      <TR><TD><p style='color:{self.color};'>{self.name}</p><TD>{self.desc}</TR>\n")
        
incomplete = Flag('MAGENTA','magenta','#FF00FF','Reduction has yet to be started due to incomplete transfer')
stopped = Flag('RED','red','#FF0000','Stopped Reduction')
NoExp = Flag('MAROON','Maroon','#FFA500','Stopped Reduction for NO GOOD EXPOSURES')
Error_warn = Flag('ORANGE','DarkOrange','#FF8C00','Pipeline ran with errors/warnings')
running = Flag('YELLOW','Gold', '#FFD700','Pipeline is still running')
NoRedux = Flag('BLUE','blue','#0000FF','No Reductions')
NoIssues = Flag('GREEN','green','#008000','No issues')


class Crash_log:
    def __init__(self, step, error,msg=None,line=None, color=Error_warn):
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
                    msg='No Sky Fibers Found', color=stopped),
          Crash_log('spDiag2d','ABORT: No good flats (saturated?)', color=stopped),
          Crash_log('spDiag2d','SPCALIB: .*: .* paired with no arc', color=stopped),
          Crash_log('spDiag2d','ABORT: Reject science as too bright: 25-th-percentile =',
                    msg='Reject Bright Science', color=stopped),
          Crash_log('spDiag2d','SKYSUBTRACT:.*: Discarding .*(fractional) of the sky pixels as bad',
                    msg='Failed Sky Subtraction', line = -1, color=stopped),
          Crash_log('spDiag2d','FITSPECTRARESOL: .*: Calculating the spectra resolution',
                    msg='Failed FITSPECTRARESOL', line = 1, color=stopped),
          Crash_log('spDiag2d','EXTRACT_BUNDLE_IMAGE: .*: sigmasize:',
                    msg='Failure Extracting Exposure', line = 1, color=stopped),
          Crash_log('spDiag2d','FITMEANX: .*:',msg='Failure in Sky Line Identification',
                    line = 1, color=stopped),
          Crash_log('spDiagcomb','RM_SPFLUX_V5:.*: USING XYFIT', color=stopped,
                    msg='SpectroPhoto Calibration Failure', line = 1),
          Crash_log('spDiagcomb','RM_SPCOMBINE_V5: ABORT: No exposures with SCORE > 0',
                    msg='No Good Exposures', color=NoExp),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: Rejected  .* of  .* std stars',
                    msg='Failure Combining Exposures', line = 1, color=stopped),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: ABORT: No good fluxing stars!',
                    color=Error_warn, msg='ABORT: No good fluxing stars!'),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: WARNING: Already rejected .* of  .* std stars',
                    color=Error_warn),
          Crash_log('spDiag1d','ZCOMPUTE: .*',msg='Failure in COMPUTECHI2 for ZFIND',
                    line = 1, color=stopped),
          Crash_log('spDiag1d','ZFIND: .*',msg='Failure in COMPUTECHI2 for ZFIND',
                                  line = 1, color=stopped),
          Crash_log('run_spTrace','Killed', msg='Failed run_spTrace', color=stopped),
          Crash_log('spAll','fieldmerge: EXITING!!', color=stopped),
          Crash_log('spSpec_reformat', 'read_spAll: ERROR: Missing .*',
                    msg='Failed spSpec_reformat: missing spAll field', color=stopped)]

py_err = [Crash_log(None,'exception:',
                     msg='Failed {step}', color=stopped),
                     Crash_log('spAll','fieldmerge: No valid spAll entries', color='red')]

noerr_cal_b = Crash_log('spDiag2d','SPCALIB: b.*: .* paired with arc',
                        msg='No error', color=NoIssues)
noerr_cal_r = Crash_log('spDiag2d','SPCALIB: r.*: .* paired with arc',
                        msg='No error', color=NoIssues)

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
            return(running.color, f'<b>{key}</b> not yet run')
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
            return(NoIssues, None)
        elif key == 'run_spTrace' and '.e.log' in file:
            with open(file) as f:
                lines = f.readlines()
                lines.reverse()
                
                for i, line in enumerate(lines):
                    for err in errors:
                        msg = err.check(i,line,key)
                        if msg is not None:
                            return(stopped,msg)
            return(NoIssues,None)
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
                return(running, None)
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
            return(running, None)
    return(running,None)
class LogCheck:
    def __init__(self, topdir, run2d, run1d, field, mjd, dither='F',
                 epoch=False, custom = None, mjd1d=None,obs=None):
        self.topdir = topdir
        self.run2d = run2d
        self.run1d = run1d
        self.field = f2s(field)
        self.mjd = mjd
        self.mjd1d = mjd1d
        self.custom = custom
        self.dither = dither
        self.epoch = epoch
        self.obs = obs.lower()

    def html(self, fbase=[], exts =None):
        rs = ''
        note = []
        colors = []
        top2d  = ptt.join(self.topdir, self.run2d)
        for i, fb in enumerate(fbase):
            cs = False if self.custom is None else True
            
            fd = field_dir(top2d, self.field, custom = cs)
            ed = 'epoch' if self.epoch else ''
            file = ptt.join(fd,ed,fb.format(field=self.field, mjd=self.mjd,
                                            custom=self.custom, mjd1d=self.mjd1d,
                                            obs = self.obs))
            if ptt.splitext(file)[-1] == '.log':
                gf = f'log'
                bf = f"<s style='color:{running.color};'>log</s> "
            else:
                if ptt.splitext(file)[-1] == '.gz':
                    ext = ptt.splitext(ptt.splitext(file)[0])[-1]
                else:
                    ext = ptt.splitext(file)[-1]
                if exts is not None:
                    ext = exts[i]
                ext = ext.replace('.','')
                gf = f'{ext}'
                bf = f"<s style='color:{running.color};'>{ext}</s> "
            flag = NoIssues
            
            if ptt.splitext(file)[-1] == '.log':
                flag, tnote = parse_log(file, custom=self.custom)
                if 'v6_1' in self.run2d:
                    if 'spCalib_QA' in fb:
                        flag = NoIssues
                        bf = f"<s style='color:{NoIssues.color};'>log</s> "
                        tnote = None
                if tnote is not None:
                    note.append(tnote)
                
            if ptt.exists(file):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file)+f" style='color:{flag.color};'>"+gf+"</A> "
                    colors.append(flag)
                    continue
                else:
                    with open(file, 'rb') as ff:
                        if len(ff.readlines()) > 100:
                            rs = rs + "<A HREF="+chpc2html(file)+f" style='color:{flag.color};'>"+gf+"</A> "
                            colors.append(flag)
                            continue
            elif '*' in file:
                if len(glob.glob(file)) > 0:
                    rs = rs + "<A HREF="+chpc2html(ptt.dirname(file))+f" style='color:{flag.color};'>"+gf+"</A> "
                    colors.append(flag)
                    continue
                    
            if ptt.exists(file.replace('.pdf','.ps')):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file.replace('.pdf','.ps'))+f" style='color:{flag.color};'>"+gf+"</A> "
                    colors.append(flag)
                    continue
                elif 'spDiagcomb' in fbase[0]:
                    if colors[0] == NoIssues:
                        color = NoIssues.color
                        bf = bf.replace(f'color:{running.color}',f'color:{NoIssues.color}')
            rs = rs + bf
            colors.append(bf)
        if self.dither == 'T':
            rs = (rs.replace(f'color:{stopped.color}',f'color:{stopped.code}')
                    .replace(f'color:{running.color}',f'color:{running.code}'))
        if f'color:{stopped.color}' in rs:
            rs  = (rs.replace(f'color:{running.color}',f'color:{stopped.color}')
                     .replace(f'color:{NoIssues.color}',f'color:{stopped.color}'))
        if 'redux-' in rs:
            rs = rs.replace('<A','<A class="redux"')
            
        return(rs, note)


def CheckRedux(topdir, run2d, run1d, field, mjd, obs, dither = 'F', epoch=False,
               plan2d=None, custom = None, mjd1d = None):
    
    lc = LogCheck(topdir, run2d, run1d, field, mjd, dither = dither,
                  epoch=epoch, custom=custom, mjd1d=mjd1d, obs=obs.lower())
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
    
    cs = False if custom is None else True
    fd = field if custom is None else field
    mj = mjd   if custom is None else mjd1d
    spec_dir = field_spec_dir(topdir, run2d, fd, mj, epoch=epoch,
                              custom = cs, custom_name = custom)
    img_dir = field_png_dir(topdir,run2d,run1d,fd,mj,epoch=epoch,
                              custom = cs, custom_name = custom)
    fd = field_dir(ptt.join(topdir,run2d),fd,custom=cs)
    if epoch:
        fd = ptt.join(fd,'epoch')
    spec_dir = ptt.relpath(spec_dir, start = fd)
    img_dir = ptt.relpath(img_dir, start = fd)

    fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd}.log',
                                                      ptt.join(spec_dir,f'spAll-{field}-{mjd}.fits.gz')])
    fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd}.log',
                                                      img_dir, spec_dir],
                                                      exts=['.log','.png','.fits'])

    if not epoch and custom is None:
        fmjd['SpCalib'],    note['SpCalib'] = lc.html(['spCalib_QA-'+run2d+'-{field}-{mjd}.log',
                                                       'spCalib_QA-'+run2d+'-{field}-{mjd}.pdf'])

    fmjd['Note'] = []
    nep = False
    for key in note:
        if note[key] is not None:
            note[key] = np.unique(np.asarray(note[key])).tolist()
            if 'No Good Exposures' in ', '.join(note[key]):
                nep = True
            if nep is True:
                fmjd[key] = (fmjd[key].replace(f'color:{stopped.color}',f'color:{NoExp.color}')
                                      .replace(f'color:{running.color}',f'color:{NoExp.color}'))
            fmjd['Note'].append(', '.join(note[key]))
    fmjd['Note'] = list(dict.fromkeys(fmjd['Note'])) #list(set(fmjd['Note']))
    try:
        fmjd['Note'].remove('')
    except:
        pass
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
    
    fdr = '*' if custom is None else custom
    cs = False if custom is None else True
    if custom is not None:
        custom_base = custom
        if obs[0].lower() == 'lco':
            custom = custom+'_lco'
        elif obs[0].lower() == 'apo':
            custom = custom+'_apo'
        fdr = custom
    fd = field_dir(ptt.join(topdir, run2d),fdr, custom = cs)


    if redux is None or len(redux) == 0:
        if epoch:
            redux = glob.glob(ptt.join(fd,'epoch', f'spPlancombepoch-??????-{mjd}.par'))
        elif custom is not None:
            redux = glob.glob(ptt.join(fd, f'spPlanCustom-{custom}-{mjd}.par'))
        else:
            redux = glob.glob(ptt.join(fd, f'spPlan2d-??????-{mjd}.par'))

    if epoch:
        plans = glob.glob(ptt.join(fd,'epoch', f'spPlancombepoch-??????-{mjd}.par'))
    elif custom is not None:
        plans = glob.glob(ptt.join(fd, f'spPlanCustom-{custom}-{mjd}.par'))
    else:
        plans = glob.glob(ptt.join(fd, f'spPlan2d-??????-{mjd}.par'))


    html = pd.DataFrame()
    rlogs = []
    for r in redux:
        if epoch:
            plan = r.replace('redux-','spPlancombepoch-')
        elif custom is not None:
            plan = r.replace('redux-', 'spPlanCustom-')
        else:
            plan = r.replace('redux-','spPlan2d-')#+'.par'
        if '.par' not in plan:
            plan = plan+'.par'
        if 'redux-' not in r:
            if epoch:
                r = r.replace('spPlancombepoch-','redux-').replace('.par','')
            elif custom is not None:
                r = r.replace('spPlanCustom-','redux-').replace('.par','')
            else:
                r = r.replace('spPlan2d-','redux-').replace('.par','')
        rlogs.extend([r+'.o',r+'.e'])
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
                SOS_log = ptt.abspath(ptt.join(getenv(sos_dir, default=''),f"{mjd}",f"logfile-{mjd}.html"))
                SOS_log = f"<a HREF={chpc2html(SOS_log)}>Log</a>" if ptt.exists(SOS_log) else "N/A"
            body.append(f"{ob.upper()} SOS: {SOS_log}")

            if ptt.exists(SOS_log):
                valid = check_hash(ptt.abspath(ptt.join(getenv(sos_dir, default=''),f"{mjd}")))
                if valid:
                    body.append(f"{ob.upper()} SOS Tranfer: Complete")
                elif ptt.exists(ptt.abspath(ptt.join(getenv(sos_dir, default=''),f"{mjd}",f"{mjd}.sha1sum"))):
                    body.append(f"{ob.upper()} SOS Tranfer: Failed")
                else:
                    pass

            transferlog_json = ptt.join(getenv('DATA_ROOT', default=''),"staging/{obs}/atlogs/{mjd}/{mjd}_status.json")
            nightlogs = ptt.join(getenv('DATA_ROOT', default=''),"staging/{obs}/reports/mos/{th}")
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
            flag,_ = parse_log(spTrace)
            spTracee = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.e.log"))
            flag2,_ = parse_log(spTracee)
            if flag2 != NoIssues:
                flag = stopped
            if run2d == 'master':
                if mjd < 60418:
                    flag = NoIssues
            if epoch:
                flag = NoIssues
            spTrace = f"<a HREF={chpc2html(spTrace)} style='color:{flag.color};'>o.log</a>" if ptt.exists(spTrace) else "N/A"

            spTrace1 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.e.log"))
            spTrace1 = f"<a HREF={chpc2html(spTrace1)} style='color:{flag.color};'>e.log</a>" if ptt.exists(spTrace1) else "N/A"

            spTrace2 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"arcs_{mjd}_{ob.lower()}.html"))
            spTrace2 = f"<a HREF={chpc2html(spTrace2)} style='color:{flag.color};'>Plots</a>" if ptt.exists(spTrace2) else "N/A"
            body.append(f"{ob.upper()} spTrace: {spTrace} {spTrace1} {spTrace2}")
    
        else:
            reduxb = ptt.abspath(ptt.join(field_dir(ptt.join(topdir, run2d),
                                                    custom, custom=True),
                                          f'redux_{custom}-{mjd}'))
            flag, _ = parse_log(reduxb.replace('redux_','spDiagcomb-')+'.log',custom=custom)
            reduxo = reduxb+'.o'
            reduxo = f"<a HREF={chpc2html(reduxo)} style='color:{flag.color};'>o</a>"
            reduxe = reduxb+'.e'
            reduxe = f"<a HREF={chpc2html(reduxe)} style='color:{flag.color};'>e</a>"
            body.append(f"Redux Coadd: {reduxo} {reduxe}")
    # spAll
    if epoch:
        sd = 'epoch'
        sf = ''
    elif custom is not None:
        sd = custom_base
        sf =f'-{custom_base}'
    else:
        sd = 'daily'
        sf = ''
    
    spAll = ptt.join(topdir,run2d,'summary',f'{sd}',f'spAll-{run2d}{sf}.fits.gz')
    if ptt.exists(spAll):
        if email:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> ({time.ctime(ptt.getmtime(spAll))})"
        else:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> <span id='spall'></span>"
        body.append(spallh)
    spAll = ptt.join(topdir,run2d,'summary',f'{sd}',f'spAll-lite-{run2d}{sf}.fits.gz')
    if ptt.exists(spAll):
        if email:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> ({time.ctime(ptt.getmtime(spAll))})"
        else:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> <span id='spall-lite'></span>"
        body.append(spallh)

    # fieldlist
    if custom is None:
        flist = ptt.join(topdir,run2d,'summary',sd,f'fieldlist-{run2d}.fits')
        if ptt.exists(flist):
            if email:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (fits)</a> ({time.ctime(ptt.getmtime(flist))})"
            else:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (fits)</a> <span id='fieldlistfits'></span>"
            body.append(flisth)
        flist = ptt.join(topdir,run2d,'summary',sd,f'fieldlist.html')
        if ptt.exists(flist):
            if email:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (html)</a> ({time.ctime(ptt.getmtime(flist))})"
            else:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (html)</a> <span id='fieldlisthtml'></span>"
            body.append(flisth)
        
    body[-1] = body[-1]+"</h3>"
    body.append(html.to_html(index=False, escape=False, justify="center").replace('<td>', '<td align="center">'))

    return(body, rlogs)

def daily_log_to_file(obs, mjd, topdir=None, run2d=None, run1d=None, redux=None,
                      html_log=None, rlogs=None, summary=True, epoch=False, custom = None):
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
            body, rlogs = daily_log_html(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d,
                                        redux=redux, epoch=epoch, custom = custom)
        else:
            body = html_log

        if rlogs is not None:
            for r in rlogs:
                dir_ = ptt.join(outdir,'redux_logs',ptt.basename(ptt.dirname(r)))
                if not ptt.exists(dir_):
                    makedirs(dir_)
                try:
                    symlink(r, ptt.join(dir_, ptt.basename(r))+'.log')
                except:
                    pass

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
    cmd.append("function createAlternativeLink(primaryLink) {")
    cmd.append("    var parts = primaryLink.split('/');")
    cmd.append("    var filep = parts[parts.length-1].split('.');")
    cmd.append("    filep = [filep[filep.length - 2],filep[filep.length - 1],'log'].join('.');")
    cmd.append("    console.log(parts)")
    cmd.append("    parts = ['redux_logs',parts[parts.length - 2],filep];")
    cmd.append("    var alternativeLink = parts.join('/');")
    cmd.append("    return alternativeLink;")
    cmd.append("}")
    cmd.append("")
    cmd.append("function checkMainLink(element) {")
    cmd.append("    var mainLinkUrl = element.href;")
    cmd.append("    var altUrl = createAlternativeLink(mainLinkUrl);")
    cmd.append("    console.log(altUrl);")
    cmd.append("    // Perform an AJAX request to check if the main link is accessible")
    cmd.append("    var xhr = new XMLHttpRequest();")
    cmd.append("    xhr.open('HEAD', altUrl);")
    cmd.append("    xhr.onload = function() {")
    cmd.append("        if (xhr.status === 200) {")
    cmd.append("            console.log('Alt link is accessible.');")
    cmd.append("            console.log(altUrl);")
    cmd.append("            element.href = altUrl;")
    cmd.append("        } else {")
    cmd.append("            console.log('Alt link is not accessible. Redirecting to main link.');")
    cmd.append("        };")
    cmd.append("    };")
    cmd.append("    xhr.onerror = function() {")
    cmd.append("        console.error('Error checking alt link status.');")
    cmd.append("    };")
    cmd.append("    xhr.send();")
    cmd.append("}")
    cmd.append("")
    cmd.append("document.addEventListener('DOMContentLoaded', function() {")
    for filep in [f'spAll-{run2d}.fits.gz',f'spAll-lite-{run2d}.fits.gz',
                         f'fieldlist-{run2d}.fits',f'fieldlist.html']:
        if epoch:
            sum_dir = 'epoch'
        elif custom is not None:
            sum_dir = 'custom'
        else:
            sum_dir = 'daily'
        filep = ptt.join(topdir, run2d, 'summary',sum_dir,filep)
        if 'spAll-lite' in filep:
            filet = 'spall-lite'
        elif 'spAll' in filep:
            filet = 'spall'
        elif 'fieldlist.html' in filep:
            filet = 'fieldlisthtml'
        else:
            filet = 'fieldlistfits'
            
            
        
    cmd.append(f"    getlastmod('{chpc2html(filep)}', '{filet}');")
    cmd.append("")
    cmd.append("    var links = document.querySelectorAll('.redux');")
    cmd.append("    links.forEach(function(link) {")
    cmd.append("        checkMainLink(link);")
    cmd.append("    });")
    cmd.append("});")
    cmd.append("")
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
        f.write('   <style>BODY {font: normal 12pt Helvetica,Arial; padding-top: 80px;}</style>\n')
        f.write('<meta name="viewport" content="width=device-width, initial-scale=1.0">\n')
        f.write('   <link rel="shortcut icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" type="image/x-icon">\n')
        f.write(' </head>\n')
        f.write(' <body>\n')
        f.write('  <div style="display: flex; justify-content: space-between; align-items: center; position: fixed; '+
                              'width: 100%; top: 0; padding:5px 20px; background-color: #f1f1f1; '+
                              'box-shadow: 0px 2px 5px rgba(0, 0, 0, 0.1); z-index: 1000; box-sizing: border-box;">\n')

        f.write('   <h2><img src="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" style="height: 78; width: auto;">\n')
        if epoch:
            f.write(f'     Epoch BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')
        elif custom is not None:
            f.write(f'     {custom} BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')
        else:
            f.write(f'     Daily BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')

        f.write(" </div>\n")
        f.write(' <div style="padding-top: 50px; padding-bottom: 50px">\n')

        f.write(" <div class='row' style='display: flex; margin-left:-5px;margin-right:-5px;'>\n")
        f.write("  <div class='column' style='flex: 50%; padding: 5px;'>\n")
        f.write('    <TABLE BORDER=2; style="top: 140px;margin-top: 20px;" id="mytable";>\n')
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
                    transfer = True
                    with open(ptt.join(directory,f'{mjd}-{ob}.html')) as fl:
                        for line in fl.readlines():
                            if  f'color:{stopped.color};' in line:
                                color=stopped.color
                                break
                            if f'color:{Error_warn.color};' in line:
                                color=Error_warn.color
                            if f'color:{running.color};' in line:
                                if color not in [Error_warn.color]:
                                    color=running.color
                            if f'color:{NoExp.color};' in line:
                                if color not in [Error_warn.color,running.color]:
                                    color=NoExp.color
                            if '???' in line:
                                color=incomplete.color
                            if 'SOS: N/A' in line:
                                sos = False
                            if 'SOS Tranfer: Failed' in line:
                                transfer = False
                            if f'{ob} spTrace: N/A N/A N/A' in line:
                                if 'v6_1' not in RUN2D:
                                    sptrace = False
                            if 'spfibermap' in line:
                                redux = True
                    if not redux:
                        if nextmjd[ob] <= int(mjd):
                            color=incomplete.color
                    transferflag = ptt.join(getenv('DATA_ROOT', default=''),f"staging/{ob.lower()}/atlogs/{mjd}/transfer-{mjd}.done")
                    if ptt.exists(transferflag):
                        if sos is False and sptrace is False and color !=incomplete.color:
                            color = NoRedux.color
                    elif not transfer:
                        color = incomplete.color
                    elif int(mjd) >= 60150:
                        color=incomplete.color
                    if type(color) is not str:
                        color = color.color
                    obs_str.append(f"<A HREF='{mjd}-{ob}.html' style='color:{color};'>{ob}</A>")
                else:
                    obs_str.append(f"<A HREF='{mjd}-{ob}.html' style='color:{stopped.color};'><S style='color:{stopped.color};'>{ob}</S></A>")
            f.write(f'      <TR><TD>{mjd}<TD>{obs_str[0]}<TD>{obs_str[1]}</TR>\n')
        f.write('    </TABLE><br>\n')
        f.write('  </div>\n') # end column
        f.write("  <div class='column' style='flex: 50%; padding: 5px;'>\n")
        f.write('    <TABLE BORDER=2 id="mytable2" style="margin-top: 20px; position: sticky; top: 140px;">\n')
        f.write("      <TR><TD><b>Color</b><TD><b>Meaning</b></TR>\n")
        f.write(incomplete.key())
        f.write(stopped.key())
        f.write(NoExp.key())
        f.write(Error_warn.key())
        f.write(running.key())
        f.write(NoRedux.key())
        f.write(NoIssues.key())
        f.write('    </TABLE><br>\n')
        f.write('  </div>\n') # end column
        f.write(' </div>\n') # end row
        f.write(' </div>\n')
        f.write(' </body>\n')
        f.write(" <footer style='display:flex; justify-content:space-between; align-items:center; position:fixed;"
                                +" width:100%; bottom:0; padding:10px 20px; background-color:#f1f1f1; box-sizing: border-box;'>\n")
        f.write("    <div style=' text-align: left;flex-shrink: 0;'>")
        f.write('last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo))#+'\n')
        f.write("</div>\n")
        f.write("    <div style=' text-align: right;flex-shrink: 0;'>")
        f.write("<A HREF='error.html' style='color:red;'> Errors</A> | ")
        f.write("<A HREF='summary.html' style='color:Green;'> Summary</A></div>\n")

        f.write(' </footer>\n')
        f.write('</html>\n')
        
    try:
        rename(ptt.join(directory,'index.html.tmp'), ptt.join(directory,'index.html'))
    except:
        pass
    
    _Summary(directory, RUN2D, epoch = False, custom=None, error=True)
    _Summary(directory, RUN2D, epoch = False, custom=None)



def daily_log_email(subject, attachment, logger, obs, mjd,
                    email_file = None, topdir=None, epoch=False,
                    run2d=None, run1d=None, content=None, custom=None,
                    from_domain="chpc.utah.edu",  redux = None):
  
    
    body, _ = daily_log_html(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d,
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


def _Summary(directory, RUN2D, epoch = False, custom=None, error=False,
             mjd = '?????', obs = '???', html = True):
    logs = glob.glob(ptt.join(directory,f'{mjd}-???.html'))
    logs = [ptt.basename(x).split('-')[0] for x in logs]
    logs = np.unique(np.asarray(logs)).tolist()
    _df = None
    set_obs = obs
    for mjd in sorted(logs,reverse=True):
        obs = sorted(glob.glob(ptt.join(directory,f'{mjd}-{set_obs}.html')))
        obs = [ptt.basename(x).split('-')[1].split('.')[0] for x in obs]
        obs_str = []
        for ob in obs:
            try:
                t_ = []
                with open(ptt.join(directory,f'{mjd}-{ob}.html'), "r", encoding="utf-8") as file:
                    table = BeautifulSoup(file, "html.parser")
                spTrace_test = table.find('h3')
                spTrace_log = [s for s in ' '.join([str(s) for s in spTrace_test.contents]).split('<br/>') if 'spTrace:' in s]
                if len(spTrace_log) == 0:
                    spTFlag = None
                elif 'color: red' in spTrace_log[0] or 'color:red' in spTrace_log[0]:
                    spTFlag = 'Failure in <b>spTrace</b>'
                elif 'color:Gold' in spTrace_log[0] or 'color: Gold' in spTrace_log[0]:
                    spTFlag = '<b>spTrace</b> in Progress'
                else:
                    spTFlag = None
                for row in table.find("table").find_all("tr"):
                    cols = row.find_all(["td", "th"])
                    cols = [str(col).replace('</td>','').replace('<td align="center">','').replace('<th>','').replace('</th>','') for col in cols]  # Keep HTML content as string
                    t_.append(cols)
                t_ = pd.DataFrame(t_)
                t_.columns = t_.iloc[0] # Set the first row as the header
                t_ = t_[1:].reset_index(drop=True)
                if 'Note' not in t_.columns:
                    continue
                if spTFlag is not None:
                    t_['Note'] = t_['Note'].apply(lambda x: spTFlag + ' ,' + x if x else spTFlag)
                if 'Note' not in t_.columns:
                    continue
                if error is True:
                    t_ = t_.loc[t_.Note != '']
            except Exception as e:
                print(e)
                continue
            if _df is None:
                _df = t_.copy()
            else:
                _df = pd.concat([_df, t_], ignore_index=True, axis=0)
    if _df is not None:
        _df.index = range(len(_df), 0, -1)
        try:
            _df = _df.fillna('')
        except Exception as e:
            pass

    if html:
        _summary_html(_df, directory, RUN2D, epoch = epoch, custom=custom, error=error)
    return _df


def _summary_html(_df, directory, RUN2D, epoch = False, custom=None, error=False):
    if epoch:
        coadd = 'Epoch'
    elif custom is not None:
        coadd = f'{custom.title()}'
    else:
        coadd = 'Daily'
    if error:
        title = f'{coadd} BOSS Pipeline Error Summary: RUN2D={RUN2D}'
    else:
        title = f'{coadd} BOSS Pipeline Summary: RUN2D={RUN2D}'
    body =     [ '<html>']
    body.append( ' <head>')
    body.append(f'  <title>{title}</title>')
    body.append( '   <style>BODY {font: normal 12pt Helvetica,Arial; padding-top:80px;}</style>\n')
    body.append( '   <style>.dataframe{top: 140px; margin-top:50px;}</style>\n')
    body.append( '   <meta name="viewport" content="width=device-width, initial-scale=1.0">')
    body.append( '   <link rel="shortcut icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" type="image/x-icon">')
    body.append( ' </head>\n')
    body.append( ' <body>\n')
    body.append( '  <div style="display: flex; justify-content: space-between; align-items: left; position: fixed; '+
            'width: 100%; top: 0; padding:5px 20px; background-color: #f1f1f1; '+
            'box-shadow: 0px 2px 5px rgba(0, 0, 0, 0.1); z-index: 1000; box-sizing: border-box;">\n')
    body.append(f'   <h2><img src="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" style="height: 78; width: auto;">{title}</h2>\n')
#    body.append(f'       {title}</h2>')
    body.append( '  </div>')
    body.append( '  <div style="padding-top: 50px; padding-bottom: 50px">\n')

    try:
        body1 = _df.to_html(index=True, escape=False, justify="center", border=2,
                            classes='dataframe').replace('<td>', '<td align="center">')
    except Exception as e:
        body1 = ''
    name = f'error.html' if error else f'summary.html'
    with open(ptt.join(directory,name), 'w') as f:
        f.write("\n".join(body))
        f.write("<br>\n".join([body1]))
        f.write(" </div>")
        f.write(" <footer style='display:flex; justify-content:space-between; align-items:center; position:fixed;"
                +" width:100%; bottom:0; padding:10px 20px; background-color:#f1f1f1; box-sizing: border-box;'>\n")
        f.write("    <div style=' text-align: left;flex-shrink: 0;'>")
        f.write('last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo))#+'\n')
        f.write("</div>\n")
        f.write("    <div style=' text-align: right;flex-shrink: 0;'>")
        if error:
            f.write("<A HREF='./' style='color:Green;'> Daily Summary</A> | ")
            f.write("<A HREF='summary.html' style='color:Green;'> Summary</A></div>\n")
        else:
            f.write("<A HREF='error.html' style='color:Red;'> Errors</A> | ")
            f.write("<A HREF='./' style='color:Green;'> Daily Summary</A></div>\n")
        f.write(' </footer>\n')
        f.write('</html>\n')


