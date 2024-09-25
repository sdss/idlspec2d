from boss_drp.utils.daily_log.Flag import *
from boss_drp.utils.daily_log.summary import _Summary, trace as trace_summary
from boss_drp import daily_dir, favicon, idlspec2d_dir

import os.path as ptt
from os import rename, getenv
import glob
import numpy as np
import datetime
from collections import OrderedDict
from jinja2 import Template
import json

from pydl.pydlutils.yanny import yanny, read_table_yanny

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

def daily_log_index(directory, RUN2D, epoch = False, custom=None, flag_noSci=False, fast_mjds = None):
    if epoch:
        title = f'Epoch BOSS Pipeline Status: {RUN2D}'
    elif custom is not None:
        titme = f'{custom} BOSS Pipeline Status: {RUN2D}'
    else:
        title = f'Daily BOSS Pipeline Status: {RUN2D}'
    logs = glob.glob(ptt.join(directory,'?????-???.html'))
    logs = [ptt.basename(x).split('-')[0] for x in logs]
    logs = np.unique(np.asarray(logs)).tolist()
    mjds_status = OrderedDict()
    nextmjd ={}
    if custom is None:
        nextmjd['APO'] = get_nextmjd(RUN2D, 'APO', nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par'))
        nextmjd['LCO'] = get_nextmjd(RUN2D, 'LCO', nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par'))
    else:
        nextmjd['APO'] = jdate
        nextmjd['LCO'] = jdate
        
    name = 'flag_noSci' if flag_noSci else 'index'

    if fast_mjds is not None:
        fast_mjds = [str(x) for x in fast_mjds]
        if ptt.exists(ptt.join(directory,name+'.json')):
            try:
                with open(ptt.join(directory,name+'.json'), 'r') as json_file:
                    mjds_status = json.load(json_file, object_pairs_hook=OrderedDict)
            except:
                pass

    for mjd in sorted(logs,reverse=True):
        if fast_mjds is not None:
            if (str(mjd) not in fast_mjds):# and (mjd in mjds_status.keys()):
                continue
        obs = sorted(glob.glob(ptt.join(directory,f'{mjd}-???.html')))
        obs = [ptt.basename(x).split('-')[1].split('.')[0] for x in obs]
        
        obs_str = []
        for ob in ['APO','LCO']:
            if ob in obs:
                color='green'
                sos=True
                sptrace=True
                rsptrace=False
                redux = False
                transfer = True
                with open(ptt.join(directory,f'{mjd}-{ob}.html')) as fl:
                    for line in fl.readlines():
                        line = " ".join(line.split())
                        if  f'color:{stopped.color};' in line:
                            color=stopped.color
                            if f'{ob} spTrace:' in line and 'N/A N/A N/A' not in line:
                                rsptrace = True
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
                        if f'{ob} spTrace: N/A N/A N/A N/A' in line:
                            if 'v6_1' not in RUN2D:
                                sptrace = False
                        elif f'{ob} spTrace:' in line:
                            rsptrace = True
                        if 'spfibermap' in line:
                            redux = True
                if not redux:
                    if rsptrace is False:
                        if nextmjd[ob] <= int(mjd):
                            color=incomplete.color
                transferflag = ptt.join(getenv('DATA_ROOT', default=''),f"staging/{ob.lower()}/atlogs/{mjd}/transfer-{mjd}.done")
                if ptt.exists(transferflag):
                    if sos is False and sptrace is False and color !=incomplete.color:
                        color = NoObs.color
                    elif sptrace is False and color !=incomplete.color and flag_noSci:
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
        mjds_status[mjd] = OrderedDict(mjd=mjd,apo=obs_str[0],lco=obs_str[1])
        
    with open(ptt.join(directory,name+'.json.tmp'), 'w') as json_file:
        jfs = json.dumps(mjds_status,indent=4)
        json_file.write(jfs)
    try:
        rename(ptt.join(directory,name+'.json.tmp'), ptt.join(directory,name+'.json'))
    except:
        pass
    mjds_status = [mjds_status[mjd] for mjd in mjds_status.keys()]

    if flag_noSci:
        key = [incomplete.key(),stopped.key(),NoExp.key(),Error_warn.key(),
               running.key(),NoObs.key(),NoRedux.key(),NoIssues.key()]
        footer = ["<A HREF='./' style='color:#008000;'> Daily Summary</A>",
                  "<A HREF='trace.html' style='color:#0000FF;'> spTrace</A>",
                  "<A HREF='error.html' style='color:Red;'> Errors</A>",
                  "<A HREF='summary.html' style='color:Green;'> Summary</A>"]
    else:
        key = [incomplete.key(),stopped.key(),NoExp.key(),Error_warn.key(),
               running.key(),NoObs.key(),NoIssues.key()]
        footer = ["<A HREF='flag_noSci.html' style='color:#48D1CC;'> flag_noSci</A>",
                  "<A HREF='trace.html' style='color:#0000FF;'> spTrace</A>",
                  "<A HREF='error.html' style='color:Red;'> Errors</A>",
                  "<A HREF='summary.html' style='color:Green;'> Summary</A>"]
    lastupdate=('last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo))

    template = ptt.join(idlspec2d_dir,'templates','html','daily_index_template.html')

    with open(ptt.join(directory,name+'.html.tmp'), 'w',encoding="utf-8") as f:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            f.write(j2_template.render(title=title, favicon=favicon,
                                       daily_status=mjds_status, key=key,
                                       lastupdate = lastupdate, footer=footer))
        
    try:
        rename(ptt.join(directory,name+'.html.tmp'), ptt.join(directory,name+'.html'))
    except:
        pass
    
    if not flag_noSci:
        _Summary(directory, RUN2D, epoch = epoch, custom=custom, error=True)
        _Summary(directory, RUN2D, epoch = epoch, custom=custom)
        if epoch is False and custom is None:
            trace_summary(directory, RUN2D)
