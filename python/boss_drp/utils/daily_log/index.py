from boss_drp.utils.daily_log.Flag import (incomplete, stopped, NoExp, Error_warn,
                                           running, NoRedux, NoIssues)
from boss_drp.utils.daily_log.summary import _Summary, trace as trace_summary
from boss_drp import daily_dir

import os.path as ptt
from os import rename, getenv
import glob
import numpy as np
import datetime

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

def daily_log_index(directory, RUN2D, epoch = False, custom=None):
    with open(ptt.join(directory,'index.html.tmp'), 'w') as f:
        f.write('<html>\n')
        f.write(' <head>\n')
        if epoch:
            f.write(f'   <title>Epoch BOSS Pipeline Status: {RUN2D}</title>\n')
        elif custom is not None:
            f.write(f'   <title>{custom} BOSS Pipeline Status: {RUN2D}</title>\n')
        else:
            f.write(f'   <title>Daily BOSS Pipeline Status: {RUN2D}</title>\n')
        f.write('   <style>BODY {font: normal 12pt Helvetica,Arial; padding-top: 80px;}</style>\n')
        f.write('<meta name="viewport" content="width=device-width, initial-scale=1.0">\n')
        f.write('   <link rel="shortcut icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" type="image/x-icon">\n')
        f.write(' </head>\n')
        f.write(' <body>\n')
        f.write('  <div style="display: flex; justify-content: space-between; align-items: center; position: fixed; '+
                              'width: 100%; top: 0; padding:5px 20px; background-color: #f1f1f1; '+
                              'box-shadow: 0px 2px 5px rgba(0, 0, 0, 0.1); z-index: 1000; box-sizing: border-box;">\n')

        f.write('   <h2><img src="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" style="height: 78; width: auto;">\n')
        if epoch:
            f.write(f'     Epoch BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')
        elif custom is not None:
            f.write(f'     {custom} BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')
        else:
            f.write(f'     Daily BOSS Pipeline Status: RUN2D={RUN2D}</h2>\n')

        f.write(" </div>\n")
        f.write(' <div style="padding-top: 50px; padding-bottom: 50px">\n')

        f.write(" <div class='row' style='display: flex; margin-left:-5px;margin-right:-5px;'>\n")
        f.write("  <div class='column' style='flex: 50%; padding: 5px;'>\n")
        f.write('    <TABLE BORDER=2; style="top: 140px;margin-top: 20px;" id="mytable";>\n')
        logs = glob.glob(ptt.join(directory,'?????-???.html'))
        logs = [ptt.basename(x).split('-')[0] for x in logs]
        logs = np.unique(np.asarray(logs)).tolist()

        nextmjd ={}
        if custom is None:
            nextmjd['APO'] = get_nextmjd(RUN2D, 'APO', nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par'))
            nextmjd['LCO'] = get_nextmjd(RUN2D, 'LCO', nextmjd_file = ptt.join(daily_dir,'etc','nextmjd.par'))
        else:
            nextmjd['APO'] = jdate
            nextmjd['LCO'] = jdate

        for mjd in sorted(logs,reverse=True):
            obs = sorted(glob.glob(ptt.join(directory,f'{mjd}-???.html')))
            obs = [ptt.basename(x).split('-')[1].split('.')[0] for x in obs]
            
            obs_str = []
            for ob in ['APO','LCO']:
                if ob in obs:
                    color='green'
                    sos=True
                    sptrace=True
                    redux = False
                    transfer = True
                    with open(ptt.join(directory,f'{mjd}-{ob}.html')) as fl:
                        for line in fl.readlines():
                            if  f'color:{stopped.color};' in line:
                                color=stopped.color
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
                            if f'{ob} spTrace: N/A N/A N/A' in line:
                                if 'v6_1' not in RUN2D:
                                    sptrace = False
                            if 'spfibermap' in line:
                                redux = True
                    if not redux:
                        if nextmjd[ob] <= int(mjd):
                            color=incomplete.color
                    transferflag = ptt.join(getenv('DATA_ROOT', default=''),f"staging/{ob.lower()}/atlogs/{mjd}/transfer-{mjd}.done")
                    if ptt.exists(transferflag):
                        if sos is False and sptrace is False and color !=incomplete.color:
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
            f.write(f'      <TR><TD>{mjd}<TD>{obs_str[0]}<TD>{obs_str[1]}</TR>\n')
        f.write('    </TABLE><br>\n')
        f.write('  </div>\n') # end column
        f.write("  <div class='column' style='flex: 50%; padding: 5px;'>\n")
        f.write('    <TABLE BORDER=2 id="mytable2" style="margin-top: 20px; position: sticky; top: 140px;">\n')
        f.write("      <TR><TD><b>Color</b><TD><b>Meaning</b></TR>\n")
        f.write(incomplete.key())
        f.write(stopped.key())
        f.write(NoExp.key())
        f.write(Error_warn.key())
        f.write(running.key())
        f.write(NoRedux.key())
        f.write(NoIssues.key())
        f.write('    </TABLE><br>\n')
        f.write('  </div>\n') # end column
        f.write(' </div>\n') # end row
        f.write(' </div>\n')
        f.write(' </body>\n')
        f.write(" <footer style='display:flex; justify-content:space-between; align-items:center; position:fixed;"
                                +" width:100%; bottom:0; padding:10px 20px; background-color:#f1f1f1; box-sizing: border-box;'>\n")
        f.write("    <div style=' text-align: left;flex-shrink: 0;'>")
        f.write('last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo))#+'\n')
        f.write("</div>\n")
        f.write("    <div style=' text-align: right;flex-shrink: 0;'>")
        f.write("<A HREF='trace.html' style='color:#0000FF;'> spTrace</A> | ")
        f.write("<A HREF='error.html' style='color:#FF0000;'> Errors</A> | ")
        f.write("<A HREF='summary.html' style='color:#008000;'> Summary</A></div>\n")

        f.write(' </footer>\n')
        f.write('</html>\n')
        
    try:
        rename(ptt.join(directory,'index.html.tmp'), ptt.join(directory,'index.html'))
    except:
        pass
    
    _Summary(directory, RUN2D, epoch = False, custom=None, error=True)
    _Summary(directory, RUN2D, epoch = False, custom=None)
    trace_summary(directory, RUN2D)
