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
import smtplib
#from email.message import EmailMessage
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from pydl.pydlutils.yanny import yanny, read_table_yanny
import re
from os import getenv
import numpy as np

def chpc2html(fpath):
    return(fpath.replace("/uufs/chpc.utah.edu/common/home/sdss50","https://data.sdss5.org/sas/"))


def parse_log(file):
    complete= {'spfibermap':'Successful completion of readfibermaps',
               'spDiag2d':'Successful completion of SPREDUCE2D',
               'spDiagcomb':'Successful completion of SPCOMBINE',
               'spDiag1d':'Successful completion of SPREDUCE1D',
               'spXCSAO':'CPU time to compute RVs',
               'fieldlist':'Successful completion of fieldlist',
               'spAll':'Successful completion of fieldmerge',
               'spSpec_reformat':'Successful completion of spSpec_reformat',
               'spCalib_QA':'SpectroPhoto QA Complete',
               'run_spTrace':'Successful completion of boss_arcs_to_trace'}
    for key in complete:
        if key not in file:
            continue
        line = -2 if key == 'spfibermap' else -1
        with open(file) as f:
            last_line = f.readlines()[line]
        if complete[key] in last_line:
            if key == 'spDiagcomb':
                with open(file) as f:
                    last_line = f.readlines()[-4]
                if 'ABORT: No exposures with SCORE > 0' in last_line:
                    return('orange','No Good Exposures')
            if key == 'spDiag2d':
                teststr = re.compile('LOCATESKYLINES:.*WARNING: Maximum sky-line shift is.*(DISABLING)')
                with open(file) as f:
                    for line in f.readlines():
                        if teststr.match(line):
                            return('orange', teststr)
            return('green', None)
        else:
            return('orange', None)


class LogCheck:
    def __init__(self, topdir, run2d, run1d, field, mjd):
        self.topdir = topdir
        self.run2d = run2d
        self.run1d = run1d
        self.field = f2s(field)
        self.mjd = mjd

    def html(self, fbase=[], exts =None):
        rs = ''
        note = []
        for i, fb in enumerate(fbase):
            file = ptt.join(self.topdir, self.run2d, self.field,
                            fb.format(field=self.field, mjd=self.mjd))
            if ptt.splitext(file)[-1] == '.log':
                gf = f'log'
                bf = f"<s style='color:red;'>log</s> "
            else:
                if ptt.splitext(file)[-1] == '.gz':
                    ext = ptt.splitext(ptt.splitext(file)[0])[-1]
                else:
                    ext = ptt.splitext(file)[-1]
                if exts is not None:
                    ext = exts[i]
                ext = ext.replace('.','')
                gf = f'{ext}'
                bf = f"<s style='color:red;'>{ext}</s> "
            color = 'green'
            
            if ptt.splitext(file)[-1] == '.log':
                color, tnote = parse_log(file)
                if tnote is not None:
                    note.append(tnote)
            if ptt.exists(file):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file)+f" style='color: {color};'>"+gf+"</A> "
                    continue
            if ptt.exists(file.replace('.pdf','.ps')):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file.replace('.pdf','.ps'))+f" style='color: {color};'>"+gf+"</A> "
                    continue
            rs = rs + bf
        return(rs, note)


#    def txt(self, fbase=[], exts =None):
#        rs = ''
#        for fb in fbase:
#            file = ptt.join(self.topdir, self.run2d, self.field,
#                            fb.format(field=self.field, mjd=self.mjd))
#            if ptt.exists(file):
#                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
#                    rs = rs + "+ "
#                    continue
#            if ptt.exists(file.replace('.pdf','.ps')):
#                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
#                    rs = rs + "+ "
#                    continue
#            rs = rs + "- "
#        return(rs)

def CheckRedux(topdir, run2d, run1d, field, mjd, obs):
    lc = LogCheck(topdir, run2d, run1d, field, mjd)
    fmjd = pd.Series({}, dtype=object)
    note = {}
    fmjd['Field'] = field
    fmjd['MJD']   = mjd
    fmjd['OBS']   = obs.upper()
    fmjd['spfibermap'], note['spfibermap'] = lc.html(['spfibermap-{field}-{mjd}.log'])
    fmjd['spreduce2D'], note['spreduce2D'] = lc.html(['spDiag2d-{field}-{mjd}.log',
                                                      'spDiag2d-{field}-{mjd}.pdf'])
    fmjd['specombine'], note['specombine'] = lc.html(['spDiagcomb-{field}-{mjd}.log',
                                                      'spDiagcomb-{field}-{mjd}.pdf',
                                                      'spFluxdistort-{field}-{mjd}.pdf',
                                                      'spSN2d-{field}-{mjd}.pdf'])
    fmjd['spreduce1d'], note['spreduce1d'] = lc.html([run1d+'/spDiag1d-{field}-{mjd}.log'])
    fmjd['spXCSAO'],    note['spXCSAO']    = lc.html([run1d+'/spXCSAO-{field}-{mjd}.log'])
    fmjd['Fieldlist'],  note['Fieldlist']  = lc.html(['fieldlist-{field}-{mjd}.log'])
    fmjd['Fieldmerge'], note['Fieldmerge'] = lc.html(['spAll-{field}-{mjd}.log',
                                        '../spectra/full/{field}/{mjd}/spAll-{field}-{mjd}.fits.gz'])
    fmjd['Reformat'],   note['Reformat']   = lc.html(['spSpec_reformat-{field}-{mjd}.log',
                                                     f'../../images/{run2d}/{run1d}/'+'{field}-{mjd}'],
                                                     exts=['.log','.png'])
    
    fmjd['SpCalib'],    note['SpCalib'] = lc.html(['spCalib_QA-'+run2d+'-{field}-{mjd}.log',
                                                   'spCalib_QA-'+run2d+'-{field}-{mjd}.pdf'])

    fmjd['Note'] = []
    for key in note:
        if note[key] is not None:
            fmjd['Note'].append(','.join(note[key]))
    fmjd['Note'] = list(set(fmjd['Note']))
    try:
        fmjd['Note'].remove('')
    except:
        pass
    fmjd['Note'] = ','.join(fmjd['Note'])

    tfmjd = pd.Series({},dtype=object)
#    tfmjd['Field']       = field
#    tfmjd['MJD']         = mjd
#    tfmjd['spfibermap']  = lc.txt(['spfibermap-{field}-{mjd}.log'])
#    tfmjd['spreduce2D']  = lc.txt(['spDiag2d-{field}-{mjd}.log','spDiag2d-{field}-{mjd}.pdf'])
#    tfmjd['specombine']  = lc.txt(['spDiagcomb-{field}-{mjd}.log',
#                                   'spDiagcomb-{field}-{mjd}.pdf',
#                                   'spFluxdistort-{field}-{mjd}.pdf',
#                                   'spSN2d-{field}-{mjd}.pdf'])
#    tfmjd['spreduce1d']  = lc.txt([run2d+'/spDiag1d-{field}-{mjd}.log'])
#    tfmjd['spXCSAO']     = lc.txt([run2d+'/spXCSAO-{field}-{mjd}.log'])
#    tfmjd['Fieldlist']   = lc.txt(['fieldlist-{field}-{mjd}.log'])
#    tfmjd['Fieldmerge']  = lc.txt(['spAll-{field}-{mjd}.log',
#                                   '../spectra/full/{field}/{mjd}/spAll-{field}-{mjd}.fits.gz'])
#    tfmjd['Reformat']    = lc.txt(['spSpec_reformat-{field}-{mjd}.log',
#                                  f'../../images/{run2d}/{run2d}/'+'{field}-{mjd}'],
#                                  exts=['.log','.png'])
#    
#    tfmjd['SpCalib']     = lc.txt(['spCalib_QA-'+run2d+'-{field}-{mjd}.log',
#                                   'spCalib_QA-'+run2d+'-{field}-{mjd}.pdf'])
    

    return(fmjd, tfmjd)

def daily_log_email(subject, email_file, attachment, logger,
                    topdir, run2d, run1d, obs, mjd, content=None,
                    from_domain="chpc.utah.edu",  redux = None):
  
    obs = np.atleast_1d(obs).tolist()

    if redux is None:
        redux = glob.glob(ptt.join(topdir, run2d, '??????', f'redux-??????-{mjd}'))

    html = pd.DataFrame()
    txt = pd.DataFrame()
    for r in redux:
        plan = r.replace('redux-','spPlan2d-')+'.par'
        yplan = yanny(plan)
        hdr = yplan.new_dict_from_pairs()
        if 'OBS' in hdr.keys():
            thisobs = hdr['OBS'].lower()
        else:
            thisobs = 'apo'
        if thisobs not in obs:
            continue
        r = ptt.basename(r).split('-')
        field = r[1]
        mj = r[2]
        if int(mj) != int(mjd):
            continue
        fhtml, ftxt = CheckRedux(topdir, run2d, run1d, field, mj, thisobs)
        
        html = pd.concat([html, pd.DataFrame([fhtml])])
        txt = pd.concat([txt, pd.DataFrame([ftxt])])


    body = [f"<h3>RUN2D: {run2d}",f"Observatory: {','.join(obs)}"]
    for ob in obs:
        SOS_log = ptt.abspath(ptt.join(topdir,'..','sos',ob.lower(),f"{mjd}",f"logfile-{mjd}.html"))
        SOS_log = f"<a HREF={chpc2html(SOS_log)}>Log</a>" if ptt.exists(SOS_log) else "N/A"
        body.append(f"{ob.upper()} SOS: {SOS_log}")

        spTrace = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.o.log"))
        color,_ = parse_log(spTrace)
        spTrace = f"<a HREF={chpc2html(spTrace)} style='color: {color};'>o.log</a>" if ptt.exists(spTrace) else "N/A"

        spTrace1 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.e.log"))
        spTrace1 = f"<a HREF={chpc2html(spTrace1)} style='color: {color};'>e.log</a>" if ptt.exists(spTrace1) else "N/A"

        spTrace2 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"arcs_{mjd}_{ob.lower()}.html"))
        spTrace2 = f"<a HREF={chpc2html(spTrace2)} style='color: {color};'>Plots</a>" if ptt.exists(spTrace2) else "N/A"
        body.append(f"{ob.upper()} spTrace: {spTrace} {spTrace1} {spTrace2}")

    # spAll
    spAll = ptt.join(topdir,run2d,f'spAll-{run2d}.fits.gz')
    if ptt.exists(spAll):
        spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> ({ptt.getmtime(spAll)})"
        body.append(fspallh)
    spAll = ptt.join(topdir,run2d,f'spAll-lite-{run2d}.fits.gz')
    if ptt.exists(spAll):
        spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> ({ptt.getmtime(spAll)})"
        body.append(fspallh)

    # fieldlist
    flist = ptt.join(topdir,run2d,f'spAll-lite-{run2d}.fits.gz')
    if ptt.exists(flist):
        flisth = f"<a HREF={chpc2html(flist)}> FieldList</a> ({ptt.getmtime(flist)})"
        body.append(flisth)
        
        
    body[-1] = body[-1]+"</h3>"
    body.append(html.to_html(index=False, escape=False, justify="center").replace('<td>', '<td align="center">'))
    
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
                    #txt.to_string(index=False), "plain")
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

if __name__ == '__main__':
    subject = "test"
    topdir = '/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/redux/'
    run2d = 'master'
    run1d = 'master'
    obs = ['apo','lco']
    mjd = 60409
    redux = None
    attachment = None
    logger = None
    daily_log_email(subject, ptt.join(getenv('HOME'), 'daily', 'etc','emails'),
                    attachment, logger, topdir, run2d, run1d, obs, mjd, content=None,
                        from_domain="chpc.utah.edu",  redux = redux)
