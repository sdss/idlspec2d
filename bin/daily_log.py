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
    
jdate = int(float(astropy.time.Time(datetime.datetime.utcnow()).jd) - 2400000.5)

    
try:
    daily_dir = getenv('DAILY_DIR')
except:
    daily_dir = ptt.join(getenv('HOME'),'daily')

def chpc2html(fpath):
    return(fpath.replace("/uufs/chpc.utah.edu/common/home/sdss50","https://data.sdss5.org/sas/"))


class Crash_log:
    def __init__(self, step, error,msg=None,line=None, color='DarkOrange'):
        self.step  = step
        self.error = re.compile('.*'+error+'.*')
        self.msg = msg
        self.line = line
        self.color = color
    def check(self,i, line, step):
        if self.step != step:
            return
        if self.line is not None:
            if i > self.line:
                return
        if self.error.match(line):
            if self.msg is None:
                return(line.replace('\n',''))
            else:
                return(self.msg)
   
   
errors = [Crash_log('spDiag2d','LOCATESKYLINES:.*WARNING: Maximum sky-line shift is.*(DISABLING)'),
          Crash_log('spDiag2d','SKYSUBTRACT:.*: Discarding .*(fractional) of the sky pixels as bad',
                    msg='Failed Sky Subtraction', line = -1, color='red'),
          Crash_log('spDiag2d','ABORT: Only            0 sky fibers found',
                    msg='No Sky Fibers Found', color='red'),
          Crash_log('spDiag2d','ABORT: No good flats (saturated?)', color='red'),
          Crash_log('spDiag2d','SPCALIB: .*: .* paired with no arc', color='red'),
          Crash_log('spDiag2d','ABORT: Reject science as too bright: 25-th-percentile =',
                    msg='Reject Bright Science', color='red'),
          Crash_log('spDiag2d','FITSPECTRARESOL: .*: Calculating the spectra resolution',
                    msg='Failed FITSPECTRARESOL', line = -1, color='red'),
          Crash_log('spDiagcomb','RM_SPFLUX_V5:.*: USING XYFIT', color='red',
                    msg='SpectroPhoto Calibration Failure', line = -1),
          Crash_log('spDiagcomb','RM_SPCOMBINE_V5: ABORT: No exposures with SCORE > 0',
                    msg='No Good Exposures', color='red'),
          Crash_log('run_spTrace','Killed', msg='Failed run_spTrace', color='red')]

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
        if not ptt.exists(file):
            return('Gold', f'<b>{key}</b> not yet run')
        with open(file) as f:
            last_line = f.readlines()[line]
        if complete[key] in last_line:
            with open(file) as f:
                lines = f.readlines()
                lines.reverse()
                
                for i, line in enumerate(lines):
                    for err in errors:
                        msg = err.check(i,line,key)
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
            return('Gold', None)

class LogCheck:
    def __init__(self, topdir, run2d, run1d, field, mjd, dither='F'):
        self.topdir = topdir
        self.run2d = run2d
        self.run1d = run1d
        self.field = f2s(field)
        self.mjd = mjd
        self.dither = dither

    def html(self, fbase=[], exts =None):
        rs = ''
        note = []
        for i, fb in enumerate(fbase):
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
                color, tnote = parse_log(file)
                if tnote is not None:
                    note.append(tnote)
                
            if ptt.exists(file):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file)+f" style='color:{color};'>"+gf+"</A> "
                    continue
            if ptt.exists(file.replace('.pdf','.ps')):
                if ptt.getsize(file.replace('.pdf','.ps')) > 0:
                    rs = rs + "<A HREF="+chpc2html(file.replace('.pdf','.ps'))+f" style='color:{color};'>"+gf+"</A> "
                    continue
            rs = rs + bf
        if self.dither == 'T':
            rs = rs.replace('color:red','color:#FF0000')
            rs = rs.replace('color:Gold','color:#FFD700')
        if 'color:red' in rs:
            rs  = rs.replace('color:Gold','color:red').replace('color:green','color:red')
        return(rs, note)


def CheckRedux(topdir, run2d, run1d, field, mjd, obs, dither = 'F'):
    lc = LogCheck(topdir, run2d, run1d, field, mjd, dither = dither)
    fmjd = pd.Series({}, dtype=object)
    note = {}
    fmjd['Field']  = field
    fmjd['MJD']    = mjd
    fmjd['OBS']    = obs.upper()
    fmjd['Dither'] = dither
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
            fmjd['Note'].append(', '.join(note[key]))
    fmjd['Note'] = list(dict.fromkeys(fmjd['Note'])) #list(set(fmjd['Note']))
    try:
        fmjd['Note'].remove('')
    except:
        pass
    fmjd['Note'] = ', '.join(fmjd['Note'])
    
    return(fmjd)


def daily_log_html(obs, mjd, topdir=None, run2d=None, run1d=None, redux=None ):
    obs = np.atleast_1d(obs).tolist()
    if topdir is None:
        topdir = getenv('BOSS_SPECTRO_REDUX')
    if run2d is None:
        run2d = getenv('RUN2D')
    if run1d is None:
        run1d = getenv('RUN1D')
    if redux is None:
        redux = glob.glob(ptt.join(topdir, run2d, '??????', f'redux-??????-{mjd}'))

    plans = glob.glob(ptt.join(topdir, run2d, '??????', f'spPlan2d-??????-{mjd}.par'))

    html = pd.DataFrame()
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
        if 'DITHER' in hdr.keys():
            dither = hdr['DITHER']
        else:
            dither = '?'
        r = ptt.basename(r).split('-')
        field = r[1]
        mj = r[2]
        if int(mj) != int(mjd):
            continue
        fhtml = CheckRedux(topdir, run2d, run1d, field, mj, thisobs, dither=dither)
        
        html = pd.concat([html, pd.DataFrame([fhtml])])


    body = [f"<h3>RUN2D: {run2d}",f"Observatory: {','.join(obs)}", f"MJD: {mjd}"]
    for ob in obs:
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
        spTrace = f"<a HREF={chpc2html(spTrace)} style='color: {color};'>o.log</a>" if ptt.exists(spTrace) else "N/A"

        spTrace1 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"run_spTrace_{mjd}_{ob.upper()}.e.log"))
        spTrace1 = f"<a HREF={chpc2html(spTrace1)} style='color: {color};'>e.log</a>" if ptt.exists(spTrace1) else "N/A"

        spTrace2 = ptt.abspath(ptt.join(topdir,run2d,'trace',f"{mjd}",f"arcs_{mjd}_{ob.lower()}.html"))
        spTrace2 = f"<a HREF={chpc2html(spTrace2)} style='color: {color};'>Plots</a>" if ptt.exists(spTrace2) else "N/A"
        body.append(f"{ob.upper()} spTrace: {spTrace} {spTrace1} {spTrace2}")

    # spAll
    spAll = ptt.join(topdir,run2d,f'spAll-{run2d}.fits.gz')
    if ptt.exists(spAll):
        spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> ({time.ctime(ptt.getmtime(spAll))})"
        body.append(spallh)
    spAll = ptt.join(topdir,run2d,f'spAll-lite-{run2d}.fits.gz')
    if ptt.exists(spAll):
        spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> ({time.ctime(ptt.getmtime(spAll))})"
        body.append(spallh)

    # fieldlist
    flist = ptt.join(topdir,run2d,f'fieldlist-{run2d}.fits')
    if ptt.exists(flist):
        flisth = f"<a HREF={chpc2html(flist)}> FieldList (fits)</a> ({time.ctime(ptt.getmtime(flist))})"
        body.append(flisth)
 
    flist = ptt.join(topdir,run2d,f'fieldlist.html')
    if ptt.exists(flist):
        flisth = f"<a HREF={chpc2html(flist)}> FieldList (html)</a> ({time.ctime(ptt.getmtime(flist))})"
        body.append(flisth)
        
    body[-1] = body[-1]+"</h3>"
    body.append(html.to_html(index=False, escape=False, justify="center").replace('<td>', '<td align="center">'))

    return(body)

def daily_log_to_file(obs, mjd, topdir=None, run2d=None, run1d=None, redux=None, html_log=None):
    obs = np.atleast_1d(obs).tolist()
    if run2d is None:
        run2d = getenv('RUN2D')
    outdir = ptt.join(daily_dir, 'logs', 'Status', run2d)
    makedirs(outdir, exist_ok=True)

    for obs in obs:
        if html_log is None:
            body = daily_log_html(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d, redux=redux)
        else:
            body = html_log
    
        with open(ptt.join(outdir,f'{mjd}-{obs.upper()}.html'), 'w') as f:
            f.write("<br>\n".join(body))

    daily_log_index(outdir, run2d)

def daily_log_index(directory, RUN2D):
    with open(ptt.join(directory,'index.html.tmp'), 'w') as f:
        f.write('<html>\n')
        f.write(' <head>\n')
        f.write(f'   <title>Daily BOSS Pipeline Status: {RUN2D}</title>\n')
        f.write('   <style>BODY {font: normal 12pt Helvetica,Arial}</style>\n')
        f.write('   <link rel="shortcut icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" type="image/x-icon">\n')
        f.write(' </head>\n')
        f.write(' <body>\n')
        f.write('   <img src="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" style="height: 78; width: auto;">\n')
        f.write(f'   <h2>Daily BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')
        f.write(' <hr>\n')
        
        f.write(" <div class='row' style='display: flex; margin-left:-5px;margin-right:-5px;'>")
        f.write("  <div class='column' style='flex: 50%; padding: 5px;'>")
        f.write('    <TABLE BORDER=2; padding=15px 10px id="mytable";>\n')
        logs = glob.glob(ptt.join(directory,'?????-???.html'))
        logs = [ptt.basename(x).split('-')[0] for x in logs]
        logs = np.unique(np.asarray(logs)).tolist()

        for mjd in sorted(logs,reverse=True):
            obs = sorted(glob.glob(ptt.join(directory,f'{mjd}-???.html')))
            obs = [ptt.basename(x).split('-')[1].split('.')[0] for x in obs]
            
            obs_str = []
            for ob in ['APO','LCO']:
                if ob in obs:
                    color='green'
                    sos=True
                    sptrace=True
                    with open(ptt.join(directory,f'{mjd}-{ob}.html')) as fl:
                        for line in fl.readlines():
                            if 'color:red;' in line:
                                color='red'
                                break
                            
                            if 'color:Gold;' in line:
                                color='Gold'
                            if 'color:DarkOrange;' in line:
                                color='DarkOrange'
                            if '???' in line:
                                color='magenta'
                            if 'SOS: N/A' in line:
                                sos = False
                            if f'{ob} spTrace: N/A N/A N/A':
                                sptrace = False
                    transferflag = f"/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/staging/{ob.lower()}/atlogs/{mjd}/transfer-{mjd}.done"
                    if ptt.exists(transferflag):
                        if sos is False and sptrace is False and color !='magenta':
                            color = 'blue'
                    else:
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
        
    rename(ptt.join(directory,'index.html.tmp'), ptt.join(directory,'index.html'))
    
    


def daily_log_email(subject, attachment, logger, obs, mjd,
                    email_file = None, topdir=None,
                    run2d=None, run1d=None, content=None,
                    from_domain="chpc.utah.edu",  redux = None):
  
    
    body = daily_log_html(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d, redux=redux)
    
    daily_log_to_file(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d, redux=redux, html_log=body)
    
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

    parser.add_argument('--topdir',   type=str, default = getenv('BOSS_SPECTRO_REDUX'),
            help='Optional override value for the environment variable $BOSS_SPECTRO_REDUX')
    parser.add_argument('--run1d',    type=str, nargs='*', default=getenv('RUN1D'),
            help='Optional override value for the enviro variable $RUN1D')
    parser.add_argument('--run2d',    type=str, nargs='*', default=getenv('RUN2D'),
            help='Optional override value for the enviro variable $RUN2D')

    parser.add_argument('--email', action='store_true', help='Send each mjd status as email')
    args = parser.parse_args()

    if args.mjd is None:
        if args.mjdend is None and args.mjdstart is None:
            args.mjd = [jdate]
        elif args.mjdend is None and args.mjdstart is not None:
            args.mjd = range(args.mjdstart, jdate+1)
        elif args.mjdstart is None and args.mjdend is not None:
            args.mjd = [args.mjdend]
        else:
            args.mjd = range(args.mjdstart, args.mjdend+1)

    
    for mjd in args.mjd:
        for obs in args.obs:
            print(mjd, obs)
            if args.email:
                daily_log_email(f'Status: {mjd} {obs}', None, None, obs, mjd, content=None,
                            email_file = ptt.join(getenv('HOME'), 'daily', 'etc','emails'),
                            topdir=args.topdir, run2d=args.run2d, run1d=args.run1d,
                            from_domain="chpc.utah.edu",  redux = None)
            else:
                daily_log_to_file(obs, mjd, topdir=args.topdir, run2d=args.run2d,
                            run1d=args.run1d, redux=None, html_log=None)
