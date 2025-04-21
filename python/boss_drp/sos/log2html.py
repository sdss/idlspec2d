from boss_drp.utils import getcard
import boss_drp
from boss_drp.field import field_to_string, config_to_string
from boss_drp.oplimits import color2hex, oplimits
from boss_drp.utils.lock import lock, unlock

from collections import OrderedDict
import os
import time
import numpy as np
from astropy.io import fits
from astropy.time import Time
from jinja2 import Template
import shutil

def format_note(note):
    #Remove the first word from each line (which is the name of the
    #  IDL proc that generated the warning or abort message)
    note = ' '.join(note.strip().split(' ')[1:])
    note = note.replace('WARNING','<B><FONT COLOR="'+color2hex('YELLOW') + '">WARNING</FONT></B>')
    note = note.replace('ABORT','<B><FONT COLOR="'+color2hex('RED') + '">ABORT</FONT></B>').lstrip()
    note = f' {note}'
    return note
  
def get_value(rows, column, flavor, CCDs, elm = None, format=None, rf = False):
    vals = []
    r_vals = []
    for ccd in CCDs:
        sr = rows[rows['CAMERA'] == ccd]
        if len(sr) == 0:
            vals.append('')
            r_vals.append(None)
            continue
        temp = sr[0][column]
        if elm is not None:
            temp = temp[elm]
        r_vals.append(temp)
        field = 'PERCENTILE98' if column == 'PERCENTILE' else column
        temp = oplimits.check(flavor, field, ccd, temp, html=True, format=format)
        if format is not None:
            if type(temp) is not str:
               temp = format.format(temp)
        if type(temp) is str:
            if 'nan' in temp: temp = temp.replace('nan', '-')
        vals.append(temp)
    if rf:
        return (vals, r_vals)
    return vals

def sn2_total(sn2s, ftype, qualities, CCDs, designMode=None, ttype='TOTALSN2'):
    totals = []
    for i, ccd in enumerate(CCDs):
        totals.append(0)
        valid = 0
        for j, sn2 in enumerate(sn2s):
            if type(sn2[i]) is str:
                continue
            if sn2[i] is None:
                continue
            if np.isnan(sn2[i]):
                valid +=1
                continue
            if qualities[j].lower() != 'excellent':
                continue
            if oplimits.check('science', ftype, ccd, sn2[i]) == '':
                valid +=1
                totals[-1] += sn2[i]
        if valid == 0:
           totals[-1] = ''
    total_str = []
    for i, t in enumerate(totals):
        if (designMode is None) or (designMode == 'unknown'):
            ftype = 'TOTAL'
        else:
            ftype = designMode.upper()
        t = oplimits.check(ftype, ttype, CCDs[i], t, html=True, format='{:7.1f}')
        if type(t) is not str:
            t = '{:7.1f}'.format(t)
        if 'nan' in t: t = t.replace('nan', '-')
        if t.strip() == '0.0': t = '-'
        total_str.append(t)
    return(total_str)

def split_list_by_missing_values(missing_reference, exps):
    exps = sorted(exps)
    split_indices = np.where(np.diff(exps) > 1)[0] + 1

    # Split the list into sublists
    sublists = np.split(exps, split_indices)

    # Convert numpy arrays to lists
    sublists = [list(sublist) for sublist in sublists]
    merged = []
    for sublist in sublists:
        if not merged:
            merged.append(sublist)
        else:
            # Check missing values between last merged and current sublist
            last_sublist = merged[-1]
            missing_values = set(range(last_sublist[-1] + 1, sublist[0]))

            # If missing values are NOT in the reference, merge them
            if not missing_values & set(missing_reference):
                merged[-1].extend(sublist)
            else:
                merged.append(sublist)

    return merged


def log2html(mjd, sosdir, logfile=None, htmlfile=None, obs = None, fps=False, sdssv_sn2=False,
             sn2_15 = False, bright = False, copydir = None, verbose = True):
    if obs is None:
        obs = os.getenv('OBSERVATORY', 'apo')
    if logfile is None:
        logfile = f'logfile-{mjd}.fits'
    if htmlfile is None:
        htmlfile = f'logfile-{mjd}.html'
    mjd = int(mjd)
    
    flags = []
    if sdssv_sn2:
        flags.append('--sdssv_sn2')
    if bright:
        sn2_15 = True
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

    lfile = os.path.join(sosdir,logfile)
    if not os.path.exists(lfile):
        print(f'No valid logfile at {lfile}')
        return
    if lock(lfile, pause=2, niter = 10):
        try:
            with fits.open(lfile, memmap=False) as hdul:
                hdr = hdul[0].header.copy()
                biasdark = hdul[1].data.copy() if hdul[1].data is not None else None
                flat     = hdul[2].data.copy() if hdul[2].data is not None else None
                arc      = hdul[3].data.copy() if hdul[3].data is not None else None
                science  = hdul[4].data.copy() if hdul[4].data is not None else None
                text     = hdul[5].data.copy() if hdul[5].data is not None else None
        finally:
            unlock(lfile)
    # Load the FITS file
        
    exts = {'bias':biasdark, 'flat':flat, 'arc':arc, 'science':science, 'text': text}
    configs = []
    if verbose:
        print('Building lists of configs')
    for ext in exts.values():
        if ext is None:
            continue
        if len(ext) == 0:
            continue
        if fps:
            filt = 'CONFIG'
            configs.extend(ext['CONFIG'].tolist())
        else:
            try:
                filt='PLATE'
                configs.extend(ext['PLATE'].tolist())
            except:
                filt='FIELD'
                configs.extend(ext['FIELD'].tolist())
    configs = sorted(set(configs))

    disk_warnings = []
    configurations = []
    if verbose:
        print('Collating Data')
    for config in configs:
        config_str = config_to_string(config)
        if verbose: print(f'{filt}: {config_str}')
        design_mode = cart = designid = fieldid = None
        config_flav = {}
        notes = []
        totals = {}
        config_exp = None
        for name, ext in exts.items():
            config_rows = []  # Populate rows for this configuration
            qualities = []
            sn2s = []
            sn2s_v2 = []
            sn2s_m15 = []
            if ext is None:
                config_flav[f'{name.lower()}_rows'] = []
                continue
            sub_ext = ext[ext[filt] == config]
            if len(sub_ext) > 0:
                if design_mode is None:
                    if fps:
                        design_mode = sub_ext[0]['DESIGNMODE']
                    else:
                        design_mode = 'Plates'
                    cart = sub_ext[0]['CARTID']
                    if fps:
                        fieldid = sub_ext[0]['FIELD']
                    else:
                        try:
                            fieldid = sub_ext[0]['PLATE']
                        except:
                            fieldid = sub_ext[0]['FIELD']
                    fieldstr = field_to_string(fieldid)
                if name == 'text':
                    for note_r in sub_ext:
                        note = note_r['TEXT']
                        exp = int(note_r['EXPNUM'])
                        if (config_exp is None) or (int(exp)< config_exp):
                            config_exp = int(exp)
                        if 'SOS disk' in note:
                            note_s = ' '.join(note.split()[-6:])
                            if note_s in disk_warnings:
                                continue
                            disk_warnings.append(note_s)
                        notes.append(format_note(note))
                else:
                    flavor = sub_ext[0]['FLAVOR'].upper()
                    if designid is None:
                        if fps:
                            designid = sub_ext[0]['DESIGNID']
                        else:
                            designid = ''
                    expids = list(sub_ext['EXPNUM'])
                    for expid in sorted(set(expids)):
                        if (config_exp is None) or (int(expid)< config_exp):
                            config_exp = int(expid)
                        rows = sub_ext[sub_ext['EXPNUM'] == expid]
                        
                        expid  = rows[0]['EXPNUM']
                        expstr =str(expid).zfill(8)
                        exptime = rows[0]['EXPTIME']
                        quality = rows[0]['QUALITY']
                        qualities.append(quality)
                        UT = Time(2400000.5+rows[0]['TAI']/(24*3600), format='jd')
                        UT = f"{UT.datetime.hour:02}:{UT.datetime.minute:02} Z"
                        tflavor = rows[0]['FLAVOR'].upper()

                        cr = OrderedDict(
                            label = f"{tflavor}-{expstr}",
                            exptime = oplimits.check(tflavor, 'EXPTIME', CCDs[0],exptime , html=True, format='{:8.1f}'),
                            temp = f"{rows[0]['AIRTEMP']:6.1f}",
                            ut = UT,
                            quality = oplimits.check(tflavor, 'QUALITY', CCDs[0], quality, html=True),
                            expnum = int(expid),
                            )
                        if flavor in ['BIAS', 'DARK']:
                            fig = f"{tflavor.lower()}Plot-{expstr}.jpeg"
                            cr['percentile_link'] = f"<A HREF='../{mjd}/{fig}'>PERCENTILE98</A>"
                            cr['PERCENTILE'] =  get_value(rows, 'PERCENTILE', tflavor, CCDs, elm = 97, format='{:4d}')
                            
                        if flavor in ['FLAT']:
                            cr['ngood_link'] = 'NGOODFIBER'
                            cr['NGOODFIBER'] = get_value(rows, 'NGOODFIBER', tflavor, CCDs, format='{:4d}')
                            cr['xmid_link'] = 'XMID'
                            cr['XMID'] = get_value(rows, 'XMID', tflavor, CCDs, format='{:7.1f}')
                            cr['xsigma_link'] = 'XSIGMA'
                            cr['XMID'] = get_value(rows, 'XSIGMA', tflavor, CCDs, format='{:5.2f}')
         
                        if flavor in ['ARC']:
                            cr['wmid_label'] = 'WAVEMID'
                            cr['WAVEMID'] = get_value(rows, 'WAVEMID', tflavor, CCDs, format='{:7.1f}')
                            cr['bestcorr_label'] = 'BESTCORR'
                            cr['BESTCORR'] = get_value(rows, 'BESTCORR', tflavor, CCDs, format='{:4.2f}')
                            cr['nlamps_link'] = 'NLAMPS'
                            cr['NLAMPS'] = get_value(rows, 'NLAMPS', tflavor, CCDs, format='{:d}')
                            cr['wsigma_link'] = 'WSIGMA'
                            cr['WSIGMA'] = get_value(rows, 'WSIGMA', tflavor, CCDs, format='{:5.2f}')

                        if flavor in ['SCIENCE', 'SMEAR']:
                            cr['sps_label'] = 'SKY/SEC'
                            cr['SPS'] = get_value(rows, 'SKYPERSEC', tflavor, CCDs, format='{:8.2f}')
                            
                            jpegfile1 =  f'snplot-{mjd}-{config_str}-{expstr}.jpeg'
                            jpegfile_15= f'snplot-sdssv15-{mjd}-{config_str}-{expstr}.jpeg'
                            jpegfile_v2= f'snplot-sdssv-{mjd}-{config_str}-{expstr}.jpeg'
                            
                            if not fps:
                                cfalt = config_to_string(config,legacy=True)
                                ajpegfile1 =  f'snplot-{mjd}-{cfalt}-{expstr}.jpeg'
                                ajpegfile_15= f'snplot-sdssv15-{mjd}-{cfalt}-{expstr}.jpeg'
                                ajpegfile_v2= f'snplot-sdssv-{mjd}-{cfalt}-{expstr}.jpeg'
                                if not os.path.exists(os.path.join(sosdir,jpegfile1)):
                                    if os.path.exists(os.path.join(sosdir,ajpegfile1)):
                                        jpegfile1 = ajpegfile1
                                if not os.path.exists(os.path.join(sosdir,jpegfile_15)):
                                    if os.path.exists(os.path.join(sosdir,ajpegfile_15)):
                                        jpegfile_15 = ajpegfile_15
                                if not os.path.exists(os.path.join(sosdir,jpegfile_v2)):
                                    if os.path.exists(os.path.join(sosdir,ajpegfile_v2)):
                                        jpegfile_v2 = ajpegfile_v2


                            cr['sn2_label'] = f"<A HREF='../{mjd}/{jpegfile1}'>(S/N)^2</A>"
                            cr['sn2'], raw =  get_value(rows, 'SN2', tflavor, CCDs, format='{:7.1f}', rf = True)
                            sn2s.append(raw)
                            if (sn2_15) and (((fieldid < 100000) and (fieldid>=16000)) or (bright)): #FPS eng/commissioning Designs or bright Sn2 mode
                                cr['sn2_15label'] = f"<A HREF='../{mjd}/{jpegfile_15}'>Mag15 (S/N)^2</A>"
                                cr['sn2_15'], raw =  get_value(rows, 'SN2_15', tflavor, CCDs, format='{:7.1f}', rf = True)
                                sn2s_m15.append(raw)
                            if sdssv_sn2:
                                cr['sn2_v2label'] = f"<A HREF='../{mjd}/{jpegfile_v2}'>v2 (S/N)^2</A>"
                                cr['sn2_v2'], raw = get_value(rows, 'SN2_V2', tflavor, CCDs, format='{:7.1f}', rf = True)
                                sn2s_v2.append(raw)
            

                        config_rows.append(cr)
            
            if (name == 'science') and (len(config_rows) > 0):
                if fps:
                    if (design_mode == 'unknown') or (design_mode.strip() == '') or (design_mode is None):
                        dmode = 'fps'
                    else:
                        dmode = design_mode.replace('_eng', '').replace('_no_apogee_skies','')
                else:
                    dmode = 'plates'
                jpegfile1 =  f'snplot-{mjd}-{config_str}.jpeg'
                jpegfile_15= f'snplot-sdssv15-{mjd}-{config_str}.jpeg'
                jpegfile_v2= f'snplot-sdssv-{mjd}-{config_str}.jpeg'
                
                if not fps:
                    cfalt = config_to_string(config,legacy=True)
                    ajpegfile1 =  f'snplot-{mjd}-{cfalt}.jpeg'
                    ajpegfile_15= f'snplot-sdssv15-{mjd}-{cfalt}.jpeg'
                    ajpegfile_v2= f'snplot-sdssv-{mjd}-{cfalt}.jpeg'
                    if not os.path.exists(os.path.join(sosdir,jpegfile1)):
                        if os.path.exists(os.path.join(sosdir,ajpegfile1)):
                            jpegfile1 = ajpegfile1
                    if not os.path.exists(os.path.join(sosdir,jpegfile_15)):
                        if os.path.exists(os.path.join(sosdir,ajpegfile_15)):
                            jpegfile_15 = ajpegfile_15
                    if not os.path.exists(os.path.join(sosdir,jpegfile_v2)):
                        if os.path.exists(os.path.join(sosdir,ajpegfile_v2)):
                            jpegfile_v2 = ajpegfile_v2

                totals['sn2_label'] = f"<A HREF='../{mjd}/{jpegfile1}'>TOTAL (S/N)^2</A>"
                totals['sn2'] = sn2_total(sn2s, 'SN2', qualities, CCDs, designMode=dmode, ttype='TOTALSN2')
                
                if (sn2_15) and ((fieldid < 100000) or (bright)):
                    totals['sn2_15label'] = f"<A HREF='../{mjd}/{jpegfile_15}'>TOTAL Mag15 (S/N)^2</A>"
                    totals['sn2_15'] = sn2_total(sn2s_m15, 'SN2', qualities, CCDs, designMode=dmode, ttype='TOTALSN2_15')

                if sdssv_sn2:
                    totals['sn2_v2label'] = f"<A HREF='../{mjd}/{jpegfile_v2}'>TOTAL v2 (S/N)^2</A>"
                    totals['sn2_v2'] = sn2_total(sn2s_v2, 'SN2', qualities, CCDs, designMode=dmode, ttype='TOTALSN2_V2')

            if name != 'text':
                config_rows = sorted(config_rows, key = lambda x: x['expnum'])
                config_flav[f'{name.lower()}_rows'] = config_rows
    
        if fps:
            config_caption = f'Configuration {config} (DesignMode: {design_mode})<br>on Cart {cart} on Field #{fieldid}'
        else:
            config_caption = f'Plate {config} on Cart {cart}'
        
        configurations.append({
            "id": config,
            "min_expid": config_exp,
            "caption": config_caption,
            "design_mode": design_mode,
            "designid": designid,
            "fieldid": fieldid,
            "cart": cart,
            "notes": notes,
            "totals": totals,
            "bias_rows":config_flav['bias_rows'],
            "flat_rows":config_flav['flat_rows'],
            "arc_rows":config_flav['arc_rows'],
            "sci_rows":config_flav['science_rows']
        })

    allexps = []
    for flavor, ext in exts.items():
        if flavor == 'text':
            continue
        if ext is None:
            continue
        allexps.extend(list(ext['EXPNUM']))
    
    fconfigurations = []
    for config in configurations:
        exps = []
        for flavor in ['bias_rows','flat_rows','arc_rows','sci_rows']:
            for cr in config[flavor]:
                exps.append(cr['expnum'])
        exps = split_list_by_missing_values(allexps,exps)
        if len(exps) == 1:
            fconfigurations.append(config)
            continue

        for i, exps_g in enumerate(exps):
            cc = config.copy()
            if fps:
                cc['caption'] = (f"Configuration {config['id']} visit {i+1} "+
                                 f"(DesignMode: {config['design_mode']})<br>"+
                                 f"on Cart {config['cart']} on Field #{config['fieldid']}")
            else:
                cc['caption'] = f"Plate {config['id']} visit {i+1} on Cart {config['cart']}"
            cc['min_expid'] = min(exps_g)
            for flavor in ['bias_rows','flat_rows','arc_rows','sci_rows']:
                crs = []
                for cr in config[flavor]:
                    if cr['expnum'] in exps_g:
                        crs.append(cr)
                cc[flavor] = crs
            cc['id'] = str(config['id'])+f'_v{i}'
            cc['notes'] = [f'<b>Note</b>: Totals and Warnings/Aborts correspond to a joined version of all visits to this '+('Configuration' if fps else 'Plate')]+cc['notes']
            fconfigurations.append(cc)

    configurations = sorted(fconfigurations, key = lambda x: x['min_expid'])
    # Data for the template
    template_data = {
        "config_label": ('Configuration' if fps else 'Plate'),
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
    if os.path.exists(os.path.join(sosdir,f'Summary_{mjd}.html')):
        template_data["sos_summary_link"] = sos_summary_link
        
    arc_shift_link = f"../{mjd}/trace/{mjd}/arcs_{mjd}_{obs.lower()}.html"
    if os.path.exists(os.path.join(sosdir,'trace',f'{mjd}')):
        template_data["arc_shift_link"] = arc_shift_link

    # Save the template in the template directory
    template = os.path.join(boss_drp.idlspec2d_dir,'templates','html','SOS_log_template.html')


    reloader = "timerID=setTimeout('location.reload(true)',60000)"
    reloader = f"<BODY ONLOAD={reloader}>"
    html_curr = 'logfile-current.html'
    if verbose:
        print('Rendering data to html template format')
        
    hfile = os.path.join(sosdir,htmlfile)
    if lock(hfile, pause = 2, niter=5):
        try:
            with open(hfile, 'w', encoding="utf-8") as f:
                with open(os.path.join(sosdir,html_curr), 'w', encoding="utf-8") as fc:
                    with open(template) as template_file:
                        j2_template = Template(template_file.read())
                        rendered = j2_template.render(template_data)
                        f.write(rendered)
                        rendered = rendered.replace('<body>', reloader)
                        rendered = rendered.replace(f'<title>{obs.upper()} BOSS Spectro',
                                                    f'<title>{obs.upper()} BOSS Spectro (Current) ')
                        rendered = rendered.replace(f'<font size="+4">BOSS Spectro ',
                                                    f'<font size="+4">BOSS Spectro (Current) ')
                        fc.write(rendered)
        finally:
            unlock(hfile)
    biasdark = flat = arc = science = text = hdr = None   


    if copydir:
        if verbose:
            print('Copying htmls to '+copydir)
        os.makedirs(copydir, exist_ok=True)
        shutil.copy2(os.path.join(sosdir,htmlfile), os.path.join(copydir,htmlfile))
        shutil.copy2(os.path.join(sosdir,html_curr), os.path.join(copydir,html_curr))

