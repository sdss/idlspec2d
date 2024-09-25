
#!/usr/bin/env python
from boss_drp import idlspec2d_dir, favicon

import argparse
import glob
import os
import os.path as ptt
import numpy as np
import datetime
from jinja2 import Template

def soshtml(mjd, obs, sosdir):
    outdir='{:s}/{:d}/trace/{:d}'.format(sosdir,mjd,mjd)
    outfile = '{:s}/arcs_{:d}_{:s}.html'.format(outdir,mjd,obs)
    template = ptt.join(idlspec2d_dir,'templates','html','arc2trace_template.html')

    try:
        now = datetime.datetime.utcnow()
    except:
        now = datetime.datetime.now(datetime.UTC)
    nowstr = now.strftime("%a %b %d %H:%M:%S %Y UTC")
    
    if obs == 'lco' :
        cams = ['b2','r2']
    else :
        cams = ['b1','r1']
    
    col1_r = glob.glob(ptt.join(outdir,'sdR-??-*.fit*.png'))
    if len(col1_r) == 0:
        return
        
    tf = glob.glob(ptt.join(outdir,'spTraceFlat-??-*.fits.gz'))
    tf = [ptt.basename(t) for t in tf]
    
    yt   = [ptt.basename(x).split('.')[0].split('-')[2] for x in col1_r]
    exps = []
    for y in yt:
        row = dict(exp=y)
        row_tit = y+'<br>'
        for cam in cams:
            row[cam[0]+'Fpng'] = f'sdR-{cam}-{y}.fit.png'
            row[cam[0]+'FCpng'] = f'sdR-{cam}-{y}.fit_compare.png'
            tt = f'spTraceTab-{cam}-{y}.fits'
            ttl = f'spTraceTab-{cam}-{y}.log'
            if ptt.exists(ptt.join(outdir,tt)):
                row[f'{cam[0]}Fit'] = f'<A HREF="{tt}">f</A>'
                row[f'{cam[0]}Log'] = f'<A HREF="{ttl}">l</A>'
            elif ptt.exists(ptt.join(outdir,ttl)):
                row[f'{cam[0]}Fig'] = f'f'
                row[f'{cam[0]}Log'] = f'<A HREF="{ttl}">l</A>'
            else:
                row[f'{cam[0]}Fig'] = f'f'
                row[f'{cam[0]}Log'] = f'l'
        exps.append(row)
    
    jinja2_opts = dict(mjd=mjd, date = nowstr, tf= tf,
                        blue=cams[0], red=cams[1],
                        favicon=favicon,
                        name = f'Arc2Trace {mjd} {obs.upper()}',
                        exps=exps, obs=obs.upper())
    
    with open(outfile, 'w', encoding="utf-8") as f:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            f.write(j2_template.render(jinja2_opts))

    print('Done Updating {:s}/arcs_{:d}_{:s}.html'.format(outdir,mjd,obs))
