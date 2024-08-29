from boss_drp.utils.daily_log.html import daily_log_html
from boss_drp.utils.daily_log.index import daily_log_index
from boss_drp.utils.chpc2html import chpc2html
from boss_drp import daily_dir, favicon, idlspec2d_dir

import numpy as np
from os import getenv, makedirs, symlink
import os.path as ptt
import datetime
from jinja2 import Template

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
    sfiles = daily_log_js(outdir, topdir, run2d, epoch = epoch, custom=custom)

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
        template = ptt.join(idlspec2d_dir,'templates','html','daily_log_template.html')

        with open(ptt.join(outdir,f'{mjd}-{obs.upper()}.html'), 'w', encoding="utf-8") as f:
            name = f'{run2d} {mjd}-{obs.upper()}'
            updated =('   last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                     str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo))
            with open(template) as template_file:
                j2_template = Template(template_file.read())
                f.write(j2_template.render(summary = sfiles, body=body,
                                            favicon=favicon, name=name,
                                            updated=updated))
    
    if summary:
        daily_log_index(outdir, run2d, epoch = epoch, custom=custom)

def daily_log_js(directory, topdir, run2d, epoch=False, custom=None):
    if topdir is None:
        topdir = getenv('BOSS_SPECTRO_REDUX')

    summary =[]
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
            
            
        summary.append(dict(path=chpc2html(filep), name=filet))
    return(summary)
