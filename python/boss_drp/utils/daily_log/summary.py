from boss_drp.utils.daily_log.Flag import (incomplete, stopped, NoExp, Error_warn,
                                           running, NoRedux, NoIssues)
import glob
import os.path as ptt
import numpy as np
from bs4 import BeautifulSoup
import pandas as pd
import datetime
from collections import OrderedDict

def _Summary(directory, RUN2D, epoch = False, custom=None, error=False,
             mjd = '?????', obs = '???', html = True):
    logs = glob.glob(ptt.join(directory,f'{mjd}-???.html'))
    logs = [ptt.basename(x).split('-')[0] for x in logs]
    logs = np.unique(np.asarray(logs)).tolist()
    _df = None
    set_obs = obs
    for mjd in sorted(logs,reverse=True):
        obs = sorted(glob.glob(ptt.join(directory,f'{mjd}-{set_obs}.html')))
        obs = [ptt.basename(x).split('-')[1].split('.')[0] for x in obs]
        obs_str = []
        for ob in obs:
            try:
                t_ = []
                with open(ptt.join(directory,f'{mjd}-{ob}.html'), "r", encoding="utf-8") as file:
                    table = BeautifulSoup(file, "html.parser")
                spTrace_test = table.find('h3')
                spTrace_log = [s for s in ' '.join([str(s) for s in spTrace_test.contents]).split('<br/>') if 'spTrace:' in s]
                if len(spTrace_log) == 0:
                    spTFlag = None
                elif f'color:{stopped.color}' in spTrace_log[0] or f'color: {stopped.color}' in spTrace_log[0]:
                    spTFlag = 'Failure in <b>spTrace</b>'
                elif f'color:{running.color}' in spTrace_log[0] or f'color: {running.color}' in spTrace_log[0]:
                    spTFlag = '<b>spTrace</b> in Progress'
                else:
                    spTFlag = None
                for row in table.find("table").find_all("tr"):
                    cols = row.find_all(["td", "th"])
                    cols = [str(col).replace('</td>','').replace('<td align="center">','').replace('<th>','').replace('</th>','') for col in cols]  # Keep HTML content as string
                    t_.append(cols)
                t_ = pd.DataFrame(t_)
                t_.columns = t_.iloc[0] # Set the first row as the header
                t_ = t_[1:].reset_index(drop=True)
                if 'Note' not in t_.columns:
                    continue
                if spTFlag is not None:
                    t_['Note'] = t_['Note'].apply(lambda x: spTFlag + ' ,' + x if x else spTFlag)
                if 'Note' not in t_.columns:
                    continue
                if error is True:
                    t_ = t_.loc[t_.Note != '']
            except Exception as e:
                print(e)
                continue
            if _df is None:
                _df = t_.copy()
            else:
                _df = pd.concat([_df, t_], ignore_index=True, axis=0)
    if _df is not None:
        _df.index = range(len(_df), 0, -1)
        try:
            _df = _df.fillna('')
        except Exception as e:
            pass

    if html:
        if error:
            title = '{coadd} BOSS Pipeline Error Summary: RUN2D={RUN2D}'
            name = 'error.html'
            footer = ["<A HREF='trace.html' style='color:#0000FF;'> spTrace</A>",
                      "<A HREF='./' style='color:#008000;'> Daily Summary</A>",
                      "<A HREF='summary.html' style='color:#008000;'> Summary</A>"]
        else:
            title = '{coadd} BOSS Pipeline Summary: RUN2D={RUN2D}'
            name = 'summary.html'
            footer = ["<A HREF='trace.html' style='color:#0000FF;'> spTrace</A>",
                      "<A HREF='error.html' style='color:#FF0000;'> Errors</A>",
                      "<A HREF='./' style='color:#008000;'> Daily Summary</A>"]
        _summary_html(_df, directory, RUN2D, title,name, footer,
                      epoch = epoch, custom=custom)
    return _df


def _summary_html(_df, directory, RUN2D, title, name, footer, epoch = False, custom=None, rindex=True):
    if epoch:
        coadd = 'Epoch'
    elif custom is not None:
        coadd = f'{custom.title()}'
    else:
        coadd = 'Daily'
    title = title.format(coadd=coadd, RUN2D=RUN2D)
#    if error:
#        title = f'{coadd} BOSS Pipeline Error Summary: RUN2D={RUN2D}'
#    else:
#        title = f'{coadd} BOSS Pipeline Summary: RUN2D={RUN2D}'
    body =     [ '<html>']
    body.append( ' <head>')
    body.append(f'  <title>{title}</title>')
    body.append( '   <style>BODY {font: normal 12pt Helvetica,Arial; padding-top:80px;}</style>\n')
    body.append( '   <style>.dataframe{top: 140px; margin-top:50px;}</style>\n')
    body.append( '   <meta name="viewport" content="width=device-width, initial-scale=1.0">')
    body.append( '   <link rel="shortcut icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" type="image/x-icon">')
    body.append( ' </head>\n')
    body.append( ' <body>\n')
    body.append( '  <div style="display: flex; justify-content: space-between; align-items: left; position: fixed; '+
            'width: 100%; top: 0; padding:5px 20px; background-color: #f1f1f1; '+
            'box-shadow: 0px 2px 5px rgba(0, 0, 0, 0.1); z-index: 1000; box-sizing: border-box;">\n')
    body.append(f'   <h2><img src="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" style="height: 78; width: auto;">{title}</h2>\n')
    body.append( '  </div>')
    body.append( '  <div style="padding-top: 50px; padding-bottom: 50px">\n')

    try:
        body1 = _df.to_html(index=rindex, escape=False, justify="center", border=2,
                            classes='dataframe').replace('<td>', '<td align="center">')
    except Exception as e:
        body1 = ''
    with open(ptt.join(directory,name), 'w') as f:
        f.write("\n".join(body))
        f.write("<br>\n".join([body1]))
        f.write(" </div>")
        f.write(" <footer style='display:flex; justify-content:space-between; align-items:center; position:fixed;"
                +" width:100%; bottom:0; padding:10px 20px; background-color:#f1f1f1; box-sizing: border-box;'>\n")
        f.write("    <div style=' text-align: left;flex-shrink: 0;'>")
        f.write('last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo))#+'\n')
        f.write("</div>\n")
        f.write("    <div style=' text-align: right;flex-shrink: 0;'>")
        f.write(" | ".join(footer)+"</div>\n")
        f.write(' </footer>\n')
        f.write('</html>\n')


def trace(directory, RUN2D, mjd = '?????', obs = '???', html = True):
    logs = glob.glob(ptt.join(directory,f'{mjd}-???.html'))
    logs = [ptt.basename(x).split('-')[0] for x in logs]
    logs = np.unique(np.asarray(logs)).tolist()
    _df = None
    set_obs = obs
    for mjd in sorted(logs,reverse=True):
        obs = sorted(glob.glob(ptt.join(directory,f'{mjd}-{set_obs}.html')))
        obs = [ptt.basename(x).split('-')[1].split('.')[0] for x in obs]
        obs_str = []
        row = OrderedDict({'MJD':mjd,'APO':'','LCO':''})

        for ob in obs:
            try:
                t_ = []
                with open(ptt.join(directory,f'{mjd}-{ob}.html'), "r", encoding="utf-8") as file:
                    table = BeautifulSoup(file, "html.parser")
                spTrace_test = table.find('h3')
                spTrace_log = [s for s in ' '.join([str(s) for s in spTrace_test.contents]).split('<br/>') if 'spTrace:' in s]
                try:
                    row[ob] = spTrace_log[0].split('spTrace:')[1].replace('N/A','-')
                except:
                    continue
            except Exception as e:
                print(e)
                continue
        if _df is None:
            _df = pd.DataFrame([row])
        else:
            _df = pd.concat([_df, pd.DataFrame([row])], ignore_index=True, axis=0)

    title = 'BOSS spTrace Summary: RUN2D={RUN2D}'
    name = 'trace.html'
    footer = ["<A HREF='./' style='color:Green;'> Daily Summary</A>",
              "<A HREF='summary.html' style='color:Green;'> Summary</A>",
              "<A HREF='error.html' style='color:Red;'> Errors</A>"]
    if html:
        _summary_html(_df, directory, RUN2D, title,name, footer, rindex=False)
    return _df
