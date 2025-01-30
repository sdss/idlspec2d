#!/usr/bin/env python3
from boss_drp import idlspec2d_dir, favicon

from glob import glob
import os.path as ptt
from os import rename, getenv, makedirs
import time
import shutil
from jinja2 import Template


def get_mjd(html):
    filename = ptt.splitext(ptt.basename(html))[0]
    int_part = filename.split('-')[1]
    return(int(int_part))

def build_combine_html(sosdir, force=False):
    directory = ptt.join(sosdir,'combined')
    try:
        makedirs(directory, exist_ok=True)
    except:
        pass
    if not force:
        if ptt.exists(ptt.join(directory,'index.html')):
            file_time = ptt.getmtime(ptt.join(directory,'index.html'))
            if time.time()-file_time < 3600*15:
                print('No update required')
                return
    if not ptt.exists(ptt.join(directory,"sdss-new-logo.png")):
        shutil.copy(ptt.join(idlspec2d_dir, 'etc',"sdss-new-logo.png"), ptt.join(directory,"sdss-new-logo.png"))
    template = ptt.join(idlspec2d_dir,'templates','html','SOS_index_template.html')
    logfiles = sorted(glob(ptt.join(directory,'logfile-?????.html')),key=get_mjd,reverse=True)
    logfiles = [ptt.basename(x) for x in logfiles]
    with open(ptt.join(directory,'index.html.tmp'), 'w', encoding="utf-8") as f:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            f.write(j2_template.render(logfiles = logfiles, favicon=favicon))
    rename(ptt.join(directory,'index.html.tmp'), ptt.join(directory,'index.html'))
