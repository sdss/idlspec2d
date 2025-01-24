from boss_drp.utils import getcard
import boss_drp
from boss_drp.field import field_to_string

from pydl.pydlutils.yanny import read_table_yanny
from collections import OrderedDict
import os
import time
from astropy.io import fits
from astropy.time import Time
from jinja2 import Template
import shutil

def color2hex(color):
    color = color.strip().upper()
    if color == 'RED':
        return '#FF0000'
    if color == 'YELLOW':
        return '#909000'
    return 'black'

def format_note(note):
    note = note.replace('WARNING','<B><FONT COLOR="'+color2hex('YELLOW') + '">WARNING</FONT></B>')
    note = note.replace('ABORT','<B><FONT COLOR="'+color2hex('RED') + '">ABORT</FONT></B>')
    note = f'{note}'
    return note
  
numlimits = read_table_yanny(os.path.join(boss_drp.idlspec2d_dir,'examples','opLimits.par'),'SPECLIMIT')
textlimits = read_table_yanny(os.path.join(boss_drp.idlspec2d_dir,'examples','opLimits.par'),'TEXTLIMIT')

def is_number(string):
    try:
        float(string)  # Try converting to a float
        return True
    except ValueError:
        return False
        
def checklimits(flavor, field, camera, value, html=False):
    markstr = ''
    c = camera[0]+'*'
    if value is None or (isinstance(value, (list, tuple)) and len(value) == 0):
        if html:
            return '<span>'
        return markstr
    
    value = str(value)
    if not is_number(value):
        fm = textlimits[(textlimits['field'] == field) &
                        ((textlimits['camera'] == camera) |
                         (textlimits['camera'] == c) |
                         (textlimits['camera'] == '*')) &
                        ((textlimits['flavor'] == flavor) |
                         (textlimits['flavor'] == '*')) &
                        (textlimits['strval'] == value.strip())]
        if len(fm) > 0:
            markstr = fm[0]['color']
            if html:
                markstr = f'<span style="color:{color2hex(markstr)};font-weight:bold;">'
        
    else:
        fm = numlimits[(numlimits['field'] == field) &
                       ((numlimits['camera'] == camera) |
                        (numlimits['camera'] == c) |
                        (numlimits['camera'] == '*')) &
                       ((numlimits['flavor'] == flavor) |
                        (numlimits['flavor'] == '*')) &
                       (numlimits['lovalue'] <= float(value)) &
                       (float(value) <= numlimits['hivalue'])]
        if len(fm) > 0:
            if len(fm) > 1:
                if len(fm[fm['color'] == 'red']) >= 1 :
                    fm = fm[fm['color'] == 'red'][0]
                else:
                    fm = fm[0]
            else:
                fm = fm[0]
            markstr = fm['color']
            if html:
                markstr = f'<span style="color:{color2hex(markstr)};font-weight:bold;">'
    return markstr

def get_value(rows, column, flavor, CCDs, elm = None, format=None, rf = False):
    vals = []
    r_vals = []
    for ccd in CCDs:
        sr = rows[rows['CAMERA'] == ccd]
        if len(sr) == 0:
            vals.append('-')
            r_vals.append(0)
            continue
        temp = sr[0][column]
        if elm is not None:
            temp = temp[elm]
        r_vals.append(temp)

        if format:
            if 'd' in format:
                temp = format.format(int(temp)).strip()
            else:
                temp = format.format(temp).strip()
        else:
            temp = str(temp)
        if temp == 'nan': temp = '-'
        temp = checklimits(flavor, column, ccd, temp, html=True)+temp+'</span>'
        vals.append(temp)
    if rf:
        return (vals, r_vals)
    return vals

def config_to_string(config):
    # If `config` is a list or array, handle each element recursively
    if isinstance(config, (list, tuple)):
        return [config_to_string(c) for c in config]

    # Ensure `config` is non-negative
    if config < 0:
        config = 0

    # Format the value
    if config < 1000000:
        return f"{int(config):06d}"  # Format as 6-digit integer with leading zeros
    else:
        return str(config)  # Convert to string for values >= 1,000,000


def sn2_total(sn2s, type, qualities, CCDs, designMode=None):
    totals = []
    for i, ccd in enumerate(CCDs):
        totals.append(0)
        for j, sn2 in enumerate(sn2s):
            if qualities[j].lower() != 'excellent':
                continue
            if checklimits('science', type, ccd, sn2[i]) == '':
                totals[-1] += sn2[i]
    total_str = []
    for i, t in enumerate(totals):
        tr = t
        t = '{:7.1f}'.format(t).strip()
        if t == 'nan':
            t = '-'
        if (designMode is None) or (designMode == 'unknown'):
            type = 'TOTAL'
        else:
            type = designMode.upper()
        t = checklimits(type, 'TOTALSN2', CCDs[i], tr, html=True)+t+'</span>'
        total_str.append(t)
    return(total_str)

def log2html(mjd, sosdir, logfile=None, htmlfile=None, obs = None, fps=False, sdssv_sn2=False,
             sn2_15 = False, bright = False, copydir = None):
    if obs is None:
        obs = os.getenv('OBSERVATORY', 'apo')
    if logfile is None:
        logfile = f'logfile-{mjd}.fits'
    if htmlfile is None:
        htmlfile = f'logfile-{mjd}.html'
    
    flags = []
    if sdssv_sn2:
        flags.append('--sdssv_sn2')
    if bright:
        flags.append('--bright')
    if sn2_15:
        flags.append('--sn2_15')
        
    if fps:
        flags.append('--fps')
        if obs.lower() == 'apo':
            CCDs = ['b1','r1']
        else:
            CCDs = ['b2','r2']
    elif mjd > 59030:
        obs = 'apo'
        CCDs = ['b1','r1']
    else:
        obs = 'apo'
        CCDs = ['b1','r1','b2','r2']

    # Load the FITS file
    with fits.open(os.path.join(sosdir,logfile)) as hdul:
        # Assuming relevant data is in the first HDU
        hdr = hdul[0].header
        biasdark = hdul[1].data
        flat     = hdul[2].data
        arc      = hdul[3].data
        science  = hdul[4].data
        text     = hdul[5].data
        
    exts = [biasdark, flat, arc, science, text]
    configs = []
    for ext in exts:
        if len(ext) == 0:
            continue
        configs.extend(ext['CONFIG'].tolist())
    configs = sorted(set(configs))

    exts = {'bias':biasdark, 'flat':flat, 'arc':arc, 'science':science, 'text': text}

    

    disk_warnings = []
    configurations = []
    for config in configs:
        config_str = config_to_string(config)
        design_mode = cart = designid = fieldid = None
        config_flav = {}
        notes = []
        totals = {}
        for name, ext in exts.items():
            config_rows = []  # Populate rows for this configuration
            qualities = []
            sn2s = []
            sn2s_v2 = []
            sn2s_m15 = []
            sub_ext = ext[ext['CONFIG'] == config]
            if len(sub_ext) > 0:
                if design_mode is None:
                    design_mode = sub_ext[0]['DESIGNMODE']
                    cart = sub_ext[0]['CARTID']
                    fieidid = sub_ext[0]['FIELD']
                    fieldstr = field_to_string(fieldid)
                if name == 'text':
                    t_notes = sub_ext['TEXT']
                    for note in t_notes:
                        if 'SOS disk' in note:
                            note_s = ' '.join(note.split()[-6:])
                            if note_s in disk_warnings:
                                continue
                            disk_warnings.append(note_s)
                        notes.append(format_note(note))
                else:
                    flavor = sub_ext[0]['FLAVOR'].upper()
                    if designid is None:
                        designid = sub_ext[0]['DESIGNID']
                    expids = list(sub_ext['EXPNUM'])
                    for expid in sorted(set(expids)):
                        rows = sub_ext[sub_ext['EXPNUM'] == expid]
                        
                        expid  = rows[0]['EXPNUM']
                        expstr =str(expid).zfill(8)
                        exptime = rows[0]['EXPTIME']
                        quality = rows[0]['QUALITY']
                        qualities.append(quality)
                        UT = Time(2400000.5+rows[0]['TAI']/(24*3600), format='jd')
                        UT = f"{UT.datetime.hour:02}:{UT.datetime.minute:02} Z"

                        cr = OrderedDict(
                            label = f"{flavor}-{expstr}",
                            exptime = checklimits(flavor, 'EXPTIME', CCDs[0],exptime , html=True)+f"{exptime:8.1f}</span>",
                            temp = f"{rows[0]['AIRTEMP']:6.1f}",
                            ut = UT,
                            quality = checklimits(flavor, 'QUALITY', CCDs[0], quality, html=True)+f"{quality}</span>"
                            )
                        if flavor in ['BIAS', 'DARK']:
                            fig = f"{flavor.lower()}Plot-{expstr}.jpeg"
                            cr['percentile_link'] = f"<A HREF='../{mjd}/{fig}'>PERCENTILE98</A>"
                            cr['PERCENTILE'] =  get_value(rows, 'PERCENTILE', flavor, CCDs, elm = 97, format='{:7.1f}')
                            
                        if flavor in ['FLAT']:
                            cr['ngood_link'] = 'NGOODFIBER'
                            cr['NGOODFIBER'] = get_value(rows, 'NGOODFIBER', flavor, CCDs, format='{:4d}')
                            cr['xmid_link'] = 'XMID'
                            cr['XMID'] = get_value(rows, 'XMID', flavor, CCDs, format='{:7.1f}')
                            cr['xsigma_link'] = 'XSIGMA'
                            cr['XMID'] = get_value(rows, 'XSIGMA', flavor, CCDs, format='{:5.2f}')
         
                        if flavor in ['ARC']:
                            cr['wmid_label'] = 'WAVEMID'
                            cr['WAVEMID'] = get_value(rows, 'WAVEMID', flavor, CCDs, format='{:7.1f}')
                            cr['bestcorr_label'] = 'BESTCORR'
                            cr['BESTCORR'] = get_value(rows, 'BESTCORR', flavor, CCDs, format='{:4.2f}')
                            cr['nlamps_link'] = 'NLAMPS'
                            cr['NLAMPS'] = get_value(rows, 'NLAMPS', flavor, CCDs, format='{:d}')
                            cr['wsigma_link'] = 'WSIGMA'
                            cr['WSIGMA'] = get_value(rows, 'WSIGMA', flavor, CCDs, format='{:5.2f}')

                        if flavor in ['SCIENCE', 'SMEAR']:
                            cr['sps_label'] = 'SKY/SEC'
                            cr['SPS'] = get_value(rows, 'SKYPERSEC', flavor, CCDs, format='{:8.2f}')
                            
                            jpegfile1 =  f'snplot-{mjd}-{config_str}-{expstr}.jpeg'
                            jpegfile_15= f'snplot-sdssv15-{mjd}-{config_str}-{expstr}.jpeg'
                            jpegfile_v2= f'snplot-sdssv-{mjd}-{config_str}-{expstr}.jpeg'
                            
                            cr['sn2_label'] = f"<A HREF='../{mjd}/{jpegfile1}'>(S/N)^2</A>"
                            cr['sn2'], raw =  get_value(rows, 'SN2', flavor, CCDs, format='{:7.1f}', rf = True)
                            sn2s.append(raw)
                            if (sn2_15) and (fps) and ((fieidid < 100000) or (bright)):
                                cr['sn2_15label'] = f"<A HREF='../{mjd}/{jpegfile_15}'>Mag15 (S/N)^2</A>"
                                cr['sn2_15'], raw =  get_value(rows, 'SN2_15', flavor, CCDs, format='{:7.1f}', rf = True)
                                sn2s_m15.append(raw)
                            if sdssv_sn2:
                                cr['sn2_v2label'] = f"<A HREF='../{mjd}/{jpegfile_v2}'>v2 (S/N)^2</A>"
                                cr['sn2_v2'], raw = get_value(rows, 'SN2_V2', flavor, CCDs, format='{:7.1f}', rf = True)
                                sn2s_v2.append(raw)

                        config_rows.append(cr)
            
            if (name == 'science') and (len(config_rows) > 0):
                jpegfile1 =  f'snplot-{mjd}-{config_str}.jpeg'
                jpegfile_15= f'snplot-sdssv15-{mjd}-{config_str}.jpeg'
                jpegfile_v2= f'snplot-sdssv-{mjd}-{config_str}.jpeg'
                totals['sn2_label'] = f"<A HREF='../{mjd}/{jpegfile1}'>TOTAL (S/N)^2</A>"
                if fps:
                    if (design_mode == 'unknown') or (design_mode.strip() == '') or (design_mode is None):
                        dmode = 'fps'
                    else:
                        dmode = design_mode.replace('_eng', '').replace('_no_apogee_skies','')
                else:
                    dmode = 'plates'
                totals['sn2'] = sn2_total(sn2s, 'SN2', qualities, CCDs, designMode=dmode)
                
                if (sn2_15) and (fps) and ((fieidid < 100000) or (bright)):
                    totals['sn2_15label'] = f"<A HREF='../{mjd}/{jpegfile_15}'>TOTAL Mag15 (S/N)^2</A>"
                    totals['sn2_15'] = sn2_total(sn2s_m15, 'SN2', qualities, CCDs, designMode=dmode)

                if sdssv_sn2:
                    totals['sn2_v2label'] = f"<A HREF='../{mjd}/{jpegfile_v2}'>TOTAL v2 (S/N)^2</A>"
                    totals['sn2_v2'] = sn2_total(sn2s_v2, 'SN2', qualities, CCDs, designMode=dmode)

            if name != 'text':
                config_flav[f'{name.lower()}_rows'] = config_rows
        
        configurations.append({
            "id": config,
            "design_mode": design_mode,
            "designid": designid,
            "fieldid": fieidid,
            "cart": cart,
            "notes": notes,
            "totals": totals,
            "bias_rows":config_flav['bias_rows'],
            "flat_rows":config_flav['flat_rows'],
            "arc_rows":config_flav['arc_rows'],
            "sci_rows":config_flav['science_rows']
        })

    # Data for the template
    template_data = {
        "title": f"{obs.upper()} BOSS Spectro MJD={mjd}",
        "yesterday_link": f"../{mjd-1}/logfile-{mjd-1}.html",
        "yesterday": mjd-1,
        "page_heading": f"BOSS Spectro MJD {mjd}",
        "tomorrow_link": f"../{mjd+1}/logfile-{mjd+1}.html",
        "tomorrow": mjd+1,
        "idlspec2d_version": getcard(hdr, 'VERS2D', default=boss_drp.__version__, noNaN=True),
        "run_version": getcard(hdr, 'RUN2D', default=boss_drp.__version__, noNaN=True),
        "last_updated": time.strftime("%a %b %d %H:%M:%S %Y UT", time.gmtime()),
        "configurations": configurations,
        "CCDs":CCDs,
        "sn2_15": (sn2_15 or bright),
        "sn2_v2": sdssv_sn2,
        "flags": ' '.join(flags)
    }

    sos_summary_link = f"../{mjd}/Summary_{mjd}.html"
    if True: #os.path.exist(os.path.join("..",f"{mjd}")):
        template_data["sos_summary_link"] = sos_summary_link
        
    arc_shift_link = f"../{mjd}/trace/{mjd}/arcs_{mjd}_{obs.lower()}.html"
    if True:# os.path.exist(os.path.join("..",f"{mjd}/trace/{mjd}")):
        template_data["arc_shift_link"] = arc_shift_link

    # Save the template in the template directory
    template = 'SOS_log_template.html' #os.path.join(idlspec2d_dir,'templates','html','SOS_log_template.html')


    reloader = "timerID=setTimeout('location.reload(true)',60000)"
    reloader = f"<BODY ONLOAD={reloader}>"
    html_curr = 'logfile-current.html'
    with open(os.path.join(sosdir,htmlfile), 'w', encoding="utf-8") as f:
        with open(os.path.join(sosdir,html_curr), 'w', encoding="utf-8") as fc:
            with open(template) as template_file:
                j2_template = Template(template_file.read())
                rendered = j2_template.render(template_data)
                f.write(rendered)
                fc.write(rendered.replace('<body>', reloader))


    if copydir:
        shutil.copy2(os.path.join(sosdir,htmlfile), os.path.join(copydir,htmlfile))
        shutil.copy2(os.path.join(sosdir,html_curr), os.path.join(copydir,html_curr))

