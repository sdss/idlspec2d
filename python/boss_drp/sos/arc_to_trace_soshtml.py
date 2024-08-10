
#!/usr/bin/env python
import argparse
import glob
import os
import os.path as ptt
import numpy as np
import datetime

def soshtml(mjd, obs, sosdir):
    grid=None
    xt = []
    if obs == 'lco' :
        cams = ['b2','r2']
    else :
        cams = ['b1','r1']

    outdir='{:s}/{:d}/trace/{:d}'.format(sosdir,mjd,mjd)

    for cam in cams:
        xt.append(cam)
        xt.append(cam)
        
    col1_r = glob.glob(ptt.join(outdir,'sdR-??-*.fit*.png'))
    if len(col1_r) == 0:
        return
    yt   = [ptt.basename(x).split('.')[0].split('-')[2] for x in col1_r]
    
    yt = np.unique(np.asarray(yt)).tolist()
    ytitle = []
    for y in yt:
        row = []
        row_tit = y+'<br>'
        for cam in cams:
            row.append(f'sdR-{cam}-{y}.fit.png')
            row.append(f'sdR-{cam}-{y}.fit_compare.png')
            tt = f'spTraceTab-{cam}-{y}.fits'
            ttl = f'spTraceTab-{cam}-{y}.log'
            if ptt.exists(ptt.join(outdir,tt)):
                row_tit = row_tit + f' ({cam}: <A HREF="{tt}">f</A>,<A HREF="{ttl}">l</A>)'+'<br>'
            elif ptt.exists(ptt.join(outdir,ttl)):
                row_tit = row_tit + f' ({cam}: f,<A HREF="{ttl}">l</A>)'+'<br>'
            else:
                row_tit = row_tit + f' ({cam}:f, l)'+'<br>'
        ytitle.append(row_tit)
        if grid is None:
            grid = np.asarray([row])
        else:
            grid = np.vstack([grid, np.asarray(row)])
    header = []
    header.append('<H2> {:s} MJD = {:d}</H2>'.format(obs.upper(),mjd))
    plan = 'spPlanTrace-{:d}_{:s}.par'.format(mjd,obs.upper())
    header.append('<H3> <A HREF="{:s}">{:s}</A>'.format(plan,plan))

    tf = glob.glob(ptt.join(outdir,'spTraceFlat-??-*.fits.gz'))
    for t in tf:
        plan = ptt.basename(t)
        header.append('<BR><A HREF="{:s}">{:s}</A>'.format(plan,plan))
    header.append('</H3>')

    try:
        now = datetime.datetime.utcnow()
    except:
        now = datetime.datetime.now(datetime.UTC)
    nowstr = now.strftime("%a %b %d %H:%M:%S %Y UTC")
    header.append('<H4> This page last updated <B> {:s}</B></H4>'.format(nowstr))
    header = ' \n'.join(header)

    ny=len(grid)
    nx=0
    for row in grid: nx=np.max([len(row),nx])

    with open('{:s}/arcs_{:d}_{:s}.html'.format(outdir,mjd,obs),'w') as f:
        f.write('<HTML>\n')
        f.write('<BODY>\n')
        if header is not None :
            f.write(header+'<p>\n')
        f.write('<TABLE BORDER=2>\n')
        if xt is not None :
            f.write('<TR>\n')
            if ytitle is not None :
                f.write('<TD>\n')
            for ix in range(nx) :
                f.write('<TD>'+xt[ix]+'\n')
        for iy in range(ny) :
            f.write('<TR>\n')
            if ytitle is not None :
                f.write('<TD>'+ytitle[iy]+'\n')
            for ix in range(nx) :
                f.write('<TD>\n')
                f.write('<A HREF="'+grid[iy][ix]+'">'+
                        '<IMG SRC="'+grid[iy][ix]+'" WIDTH=100%></A>\n')
                f.write('</TD>\n')
        f.write('</TABLE>\n')
        f.write('</BODY>\n')
        f.write('</HTML>\n')
    print('Done Updating {:s}/arcs_{:d}_{:s}.html'.format(outdir,mjd,obs))
