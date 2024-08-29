from boss_drp.utils.daily_log.Flag import *
from boss_drp.utils.daily_log.parse_log import CheckRedux, parse_log
from boss_drp.utils.hash import check_hash
from boss_drp.utils.chpc2html import chpc2html
from boss_drp.field import field_to_string as f2s
from boss_drp import idlspec2d_dir
from boss_drp.field import field_dir, field_spec_dir, field_png_dir

from pydl.pydlutils.yanny import yanny, read_table_yanny

import numpy as np
from os import getenv
import os.path as ptt
import glob
try:
    import json
    read_json = True
except:
    read_json = False
import pandas as pd
import time
from collections import OrderedDict
from jinja2 import Template

def daily_log_html(obs, mjd, topdir=None, run2d=None, run1d=None, redux=None,
                    email=False, epoch = False, custom=None):
    obs = np.atleast_1d(obs).tolist()
    obs = [x.lower() for x in obs]
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
    
    body = {}
    body['run2d'] = run2d
    body['obs'] = ','.join([x.upper() for x in obs])
    body['mjd'] = mjd
    if custom is not None:
        body['daily'] = None
    else:
        body['daily'] = OrderedDict()
        body['daily']['OBS']  = [x.upper() for x in obs]
    for ob in obs:
        ob = ob.upper()
        if custom is None:
            SOS_log = ptt.abspath(ptt.join(topdir,'..','sos',ob.lower(),f"{mjd}",f"logfile-{mjd}.html"))
            if ptt.exists(SOS_log):
                SOS_log_link = f"<a HREF={chpc2html(SOS_log)}>Log</a>" if ptt.exists(SOS_log) else "N/A"
            else:
                sos_dir = 'BOSS_SOS_N' if ob.lower() == 'apo' else 'BOSS_SOS_S'
                SOS_log = ptt.abspath(ptt.join(getenv(sos_dir, default=''),f"{mjd}",f"logfile-{mjd}.html"))
                SOS_log_link = f"<a HREF={chpc2html(SOS_log)}>Log</a>" if ptt.exists(SOS_log) else "N/A"
            if SOS_log_link is not None:
                body['daily'][f'{ob}_soslog'] = SOS_log_link

            if ptt.exists(SOS_log):
                valid = check_hash(ptt.abspath(ptt.join(getenv(sos_dir, default=''),f"{mjd}")),verbose=False)
                if valid:
                    body['daily'][f'{ob}_checksum'] = f'{ob} SOS Transfer: Complete <br>'
                elif ptt.exists(ptt.abspath(ptt.join(getenv(sos_dir, default=''),f"{mjd}",f"{mjd}.sha1sum"))):
                    body['daily'][f'{ob}_checksum'] = '{ob} SOS Transfer: Failed <br>'
                else:
                    pass#body['daily'][f'{ob}_checksum'] = None

            transferlog_json = ptt.join(getenv('DATA_ROOT', default=''),"staging/{obs}/atlogs/{mjd}/{mjd}_status.json")
            nightlogs = ptt.join(getenv('DATA_ROOT', default=''),"staging/{obs}/reports/mos/{th}")
            if read_json:
                if ptt.exists(transferlog_json.format(obs=ob.lower(), mjd=mjd)):
                    with open(transferlog_json.format(obs=ob.lower(), mjd=mjd)) as lf:
                        nightlog = json.load(lf)['logfile']
                        if nightlog is None:
                            body['daily'][f'{ob}_nightlog'] = "<i>Missing</i>"
                        else:
                            nightlog = nightlogs.format(obs=ob.lower(),th = nightlog)
                            nightlogh = f"<a HREF={chpc2html(nightlog)}>{ptt.basename(nightlog)}</a>"
                            body['daily'][f'{ob}_nightlog'] = nightlogh
                elif len(plans) > 0:
                    body['daily'][f'{ob}_nightlog'] = "<i>Missing</i>"
                else:
                    body['daily'][f'{ob}_nightlog'] = "???"

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
            body['daily'][f'{ob}_sptrace'] = dict(log=spTrace,elog=spTrace1,plots=spTrace2)
    
        else:
            reduxb = ptt.abspath(ptt.join(field_dir(ptt.join(topdir, run2d),
                                                    custom, custom=True),
                                          f'redux_{custom}-{mjd}'))
            flag, _ = parse_log(reduxb.replace('redux_','spDiagcomb-')+'.log',custom=custom)
            reduxo = reduxb+'.o'
            reduxo = f"<a HREF={chpc2html(reduxo)} style='color:{flag.color};'>o</a>"
            reduxe = reduxb+'.e'
            reduxe = f"<a HREF={chpc2html(reduxe)} style='color:{flag.color};'>e</a>"
            body['reduxlog'] = reduxo
            body['reduxelog'] = reduxe
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
    body['summary'] = []
    spAll = ptt.join(topdir,run2d,'summary',f'{sd}',f'spAll-{run2d}{sf}.fits.gz')
    if ptt.exists(spAll):
        if email:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> ({time.ctime(ptt.getmtime(spAll))})"
        else:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll</a> <span id='spall'></span>"
        body['summary'].append(spallh)
    spAll = ptt.join(topdir,run2d,'summary',f'{sd}',f'spAll-lite-{run2d}{sf}.fits.gz')
    if ptt.exists(spAll):
        if email:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> ({time.ctime(ptt.getmtime(spAll))})"
        else:
            spallh = f"<a HREF={chpc2html(spAll)}> spAll-lite</a> <span id='spall-lite'></span>"
        body['summary'].append(spallh)
    # fieldlist
    if custom is None:
        flist = ptt.join(topdir,run2d,'summary',sd,f'fieldlist-{run2d}.fits')
        if ptt.exists(flist):
            if email:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (fits)</a> ({time.ctime(ptt.getmtime(flist))})"
            else:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (fits)</a> <span id='fieldlistfits'></span>"
            body['summary'].append(flisth)
        flist = ptt.join(topdir,run2d,'summary',sd,f'fieldlist.html')
        if ptt.exists(flist):
            if email:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (html)</a> ({time.ctime(ptt.getmtime(flist))})"
            else:
                flisth = f"<a HREF={chpc2html(flist)}> FieldList (html)</a> <span id='fieldlisthtml'></span>"
            body['summary'].append(flisth)

    try:
        html = html.sort_values(by=['MJD','Field'], ascending = [False,True], key=lambda col: col.astype(int))
    except:
        pass
    body['fmjdlog'] = html.to_html(index=False, escape=False, justify="center").replace('<td>', '<td align="center">')

    template = ptt.join(idlspec2d_dir,'templates','html','daily_log_body_template.html')

    with open(template) as template_file:
        j2_template = Template(template_file.read())
        body = j2_template.render(**body)

    return(body, rlogs)
