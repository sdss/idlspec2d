from boss_drp import idlspec2d_dir, favicon
from boss_drp.utils.daily_log.Flag import *
import glob
import os.path as ptt
import numpy as np
from bs4 import BeautifulSoup
import pandas as pd
import datetime
from collections import OrderedDict
from jinja2 import Template
import json

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

    try:
        _df = _df.sort_values(by=['MJD', 'OBS', 'Field'], ascending=[False, True, True],
                              key=lambda col: col.astype(int) if col.name in ['Field','MJD'] else col)
        _df = _df.iloc[::-1]
        _df.reset_index(drop=True, inplace=True)
        _df = _df.iloc[::-1]

    except Exception as e:
        pass
    try:
        body1 = _df.to_html(index=rindex, escape=False, justify="center", border=2,
                            classes='dataframe').replace('<td>', '<td align="center">')
    except Exception as e:
        body1 = ''
 
    try:
        _df[['Field','MJD','OBS','Dither','Note']].to_json(ptt.join(directory,name.replace('.html','.json.tmp')), orient='records', indent=4)
    except:
        _df.to_json(ptt.join(directory,name.replace('.html','.json.tmp')), orient='records', indent=4)

    try:
        rename(ptt.join(directory,name.replace('.html','.json.tmp')), ptt.join(directory,name.replace('.html','.json')))
    except Exception as e:
        print(e)
        pass
 
    template = ptt.join(idlspec2d_dir,'templates','html','daily_Summary_all_template.html')
    lastupdate = ('last updated: '+datetime.datetime.ctime(datetime.datetime.now())+' '+
                str(datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo))
    with open(ptt.join(directory,name), 'w', encoding="utf-8") as f:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            f.write(j2_template.render(title=title, favicon=favicon, table= body1,
                                      lastupdate = lastupdate, footer=footer))


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
                except Exception as e:
                    print(e)
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
