from boss_drp.utils.daily_log.Flag import *
from boss_drp.utils.chpc2html import chpc2html
from boss_drp.field import field_to_string as f2s
from boss_drp.field import Field

import pandas as pd
import os.path as ptt
import re
from collections import OrderedDict
import numpy as np
import time
import glob


def parse_log(file, custom=None):
    complete.set(custom)
    for key in complete.line:
        if key not in file:
            continue
        line = -2 if key == 'spfibermap' else -1
        if not ptt.exists(file):
            return(running, f'<b>{key}</b> not yet run')
        with open(file) as f:
            try:
                last_line = f.readlines()[line]
            except:
                last_line = ''
        if complete.line[key] in last_line:
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
                                        noerr = noerr_cal_.check(li, ll, key)
                                        if noerr == 'No error':
                                            msg = None
                            if msg is not None:
                                return(err.flag,msg)
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
                                return(err.flag,msg)
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
                                    return(err.flag,msg)
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
        for i, fb in enumerate(fbase):
            cs = False if self.custom is None else True
            
            fd = Field(self.topdir, self.run2d, self.field, custom_name = self.custom, epoch=self.epoch).dir()
            file = ptt.join(fd,fb.format(field=self.field, mjd=self.mjd,
                                            custom=self.custom, mjd1d=self.mjd1d,
                                            obs = self.obs))
            file = ptt.abspath(file)
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

            tempfile = file

            if '.png' in file:
                if not ptt.exists(file):
                    tempfile = file.replace('.png','.pdf')
                if not ptt.exists(tempfile):
                    tempfile = file

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
                    try:
                        if colors[0] == NoIssues:
                            color = NoIssues.color
                            bf = bf.replace(f'color:{running.color}',f'color:{NoIssues.color}')
                    except:
                        pass
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
        if self.custom is not None and 'redux_' in rs:
            rs = rs.replace('<A','<A class="redux"')

        return(rs, note)

def CheckRedux(topdir, run2d, run1d, field, mjd, obs, dither = 'F', epoch=False,
               plan2d=None, custom = None, mjd1d = None):
    
    lc = LogCheck(topdir, run2d, run1d, field, mjd, dither = dither,
                  epoch=epoch, custom=custom, mjd1d=mjd1d, obs=obs.lower())
    fmjd = pd.Series({}, dtype=object)
    note = OrderedDict()
    fmjd['Path'] = ''
    fmjd['Field']  = field
    if custom is not None:
        fmjd['MJD']    = mjd1d
    else:
        fmjd['MJD']    = mjd
    fmjd['OBS']    = obs.upper()
    fmjd['Dither'] = dither
    if custom is None:
        fmjd['redux'], _ = lc.html(['redux-{field}-{mjd}','redux-{field}-{mjd}.e',
                                    'redux-{field}-{mjd}.o'], exts=['s','e','o'])
    elif f'color:{running.color}' in lc.html(['redux_{field}-{mjd}'])[0]:
        fmjd['redux'], _ = lc.html(['redux_{field}-{mjd}_{mjd1d}','redux_{field}-{mjd}_{mjd1d}.e',
                                    'redux_{field}-{mjd}_{mjd1d}.o'], exts=['s','e','o'])
    else:
        fmjd['redux'], _ = lc.html(['redux_{field}-{mjd}', 'redux_{field}-{mjd}.e','redux_{field}-{mjd}.o',
                                    'redux_{field}-{mjd}_{mjd1d}','redux_{field}-{mjd}_{mjd1d}.e',
                                    'redux_{field}-{mjd}_{mjd1d}.o'], exts=['s','e','o','1s','1e','1o'])

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
        fmjd['plans'], _ = lc.html(['spPlanCustom-{custom}_{obs}-{mjd}.par'])
        

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
                                                          #'spFluxdistort-{field}-{mjd}.pdf',
                    # turned off since we have turned off flux distortion correction by default
                                                          'spSN2d-{field}-{mjd}.pdf'])
        fmjd['spreduce1d'], note['spreduce1d'] = lc.html([run1d+'/spDiag1d-{field}-{mjd}.log',
                                                          run1d+'/spDiag1d-{field}-{mjd}.pdf'])
        fmjd['spXCSAO'],    note['spXCSAO']    = lc.html([run1d+'/spXCSAO-{field}-{mjd}.log'])
        fmjd['Fieldlist'],  note['Fieldlist']  = lc.html(['fieldlist-{field}-{mjd}.log'])
    else:
        if f'color:{running.color}' in lc.html(['redux_{field}-{mjd}'])[0]:
            fmjd['specombine'], note['specombine'] = lc.html(['spDiagcomb-{field}-{mjd}_{mjd1d}.log',
                                                              'spSN2d-{field}-{mjd}_{mjd1d}.pdf'])
        else:
            fmjd['specombine'], note['specombine'] = lc.html(['spDiagcomb-{field}-{mjd}.log',
                                                              'spSN2d-{field}-{mjd1d}.pdf'])

        fmjd['spreduce1d'], note['spreduce1d'] = lc.html([run1d+'/spDiag1d-{field}-{mjd1d}.log',
                                                          run1d+'/spDiag1d-{field}-{mjd1d}.pdf'])
        fmjd['spXCSAO'],    note['spXCSAO']    = lc.html([run1d+'/spXCSAO-{field}-{mjd1d}.log'])
    
    cs = False if custom is None else True
    fd = field if custom is None else field
    mj = mjd   if custom is None else mjd1d
    fc = Field(topdir, run2d, fd, epoch = epoch, custom_name=custom)

    spec_dir = fc.spec_dir(mj)
    img_dir = fc.png_dir(run1d, mj)
    fd = fc.dir()
        
    fmjd['Field_str'] = fmjd['Field']
    fmjd['Path'] = (f'<A HREF={chpc2html(fd)} style="text-decoration: none;">&#128193;</A>')
    spec_dir = ptt.relpath(spec_dir, start = fd)
    img_dir = ptt.relpath(img_dir, start = fd)



    if custom is None:
        fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd}.log',
                                                          ptt.join(spec_dir,f'spAll-{field}-{mjd}.fits.gz')])
        fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd}.log',
                                                          img_dir, spec_dir],
                                                          exts=['.log','.png','.fits'])
        fmjd['SpCalib'],    note['SpCalib'] = lc.html(['spCalib_QA-'+run2d+'-{field}-{mjd}.log',
                                                       'spCalib_QA-'+run2d+'-{field}-{mjd}.png'])
    else:
        fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd1d}.log',
                                                          ptt.join(spec_dir,f'spAll-{field}-{mjd1d}.fits.gz')])
        fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd1d}.log',
                                                          img_dir, spec_dir],
                                                          exts=['.log','.png','.fits'])
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
