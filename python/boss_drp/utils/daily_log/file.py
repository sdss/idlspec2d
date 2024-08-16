from boss_drp.utils.daily_log.html import daily_log_html
from boss_drp.utils.daily_log.index import daily_log_index
from boss_drp.utils.chpc2html import chpc2html
from boss_drp import daily_dir
import numpy as np
from os import getenv, makedirs, symlink
import os.path as ptt
import datetime

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
    
