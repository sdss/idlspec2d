#!/usr/bin/env python3
from boss_drp.utils import Splog
from boss_drp.field import *

import os.path as ptt
from os import getenv, makedirs, rename
from astropy.io import fits
from astropy.table import Table
from glob import glob
import argparse, sys
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import matplotlib as mpl
import matplotlib.image as image
import numpy as np
from pydl.pydlutils.sdss import sdss_flagname
import warnings
from PIL import Image, PngImagePlugin
import time

splog = Splog()


rc_fonts = {
    "font.size": 12,
    "font.family":'serif',
    
    'axes.titlesize': 18,
    "axes.labelsize": 20,
    'axes.linewidth': 1.5,
    'axes.titlepad':8,
    
    "legend.fontsize": 15,
    'figure.titlesize': 18,
    'mathtext.fontset':'dejavuserif',
    'mathtext.default': 'it',
    'mathtext.bf': 'serif:bold',
    'mathtext.it': 'serif:italic',
    
    
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.top': True,
    'xtick.direction': 'in',
    'xtick.major.size': 12,
    'xtick.minor.size': 5,
    "xtick.labelsize": 20,
        
    'ytick.right': True,
    'ytick.direction': 'in',
    'ytick.major.size': 12,
    'ytick.minor.size': 5,
    "ytick.labelsize": 20,

    
}
mpl.use('Agg')
plt.ioff()
mpl.rcParams.update(rc_fonts)

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
            to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return  ~np.isfinite(y), lambda z: z.nonzero()[0]


def sdss_spec_smooth(loglam, flux, vdisp):
    sflux=flux
    if vdisp > 1.0:
        nlambda= len(loglam)
        pixsize= np.abs(np.log(10.)*2.99792e+5*(loglam[nlambda-1]-loglam[0])/nlambda)
        smoothing= vdisp/pixsize    # pixels
        npix= int(4.0*np.ceil(smoothing))*2+3
        klam= np.arange(npix,dtype=float)-(npix-1.)/2.
        kernel= np.exp(-0.5*np.power((klam/smoothing),2))/np.sqrt(2.*np.pi)/smoothing
        kernel= kernel/np.sum(kernel)
        
        if(len(kernel) > len(flux)):
            sflux= flux
        else:
            sflux= np.convolve(flux,kernel,mode='same')#, preserve_nan=True, nan_treatment='interpolate')
    else:
        sflux = flux
    return(sflux)


def rebin(data,kernel_size ):
    kernel = np.ones(kernel_size) / kernel_size
    return np.convolve(data, kernel, mode='same')
    
def read_spAll(specfull_dir, field, mjd):
    spAll_file = ptt.join(specfull_dir, 'spAll-'+field+'-'+mjd+'.fits')
    if not ptt.exists(spAll_file):
        if ptt.exists(spAll_file+'.gz'):
            spAll_file = spAll_file+'.gz'
        else:
            splog.log('ERROR: Missing '+spAll_file)
            exit()
    spAll = Table(fits.getdata(spAll_file,1))
    spAll_hdr = fits.getheader(spAll_file,1)
    try:
        SDSSC2BV = fits.getheader(spAll_file,0)['SDSSC2BV']
    except:
        SDSSC2BV = ''
    return(spAll, spAll_hdr, SDSSC2BV)

def read_spZall(sp1d_dir, field, mjd):
    spZall_file = ptt.join(sp1d_dir, 'spZall-'+field+'-'+mjd+'.fits')
    if not ptt.exists(spZall_file):
        if ptt.exists(spZall_file+'.gz'):
            spZall_file = spZall_file+'.gz'
        else:
            splog.log('ERROR: Missing '+spZall_file)
            exit()
    spZall = Table(fits.getdata(spZall_file,1))
    spZall_hdr = fits.getheader(spZall_file,1)
    
    return(spZall, spZall_hdr)


def read_spZline(specfull_dir, field, mjd):
    spZline_file = ptt.join(specfull_dir, 'spAllLine-'+field+'-'+mjd+'.fits')
    #spZline_file = ptt.join(sp1d_dir, 'spZline-'+field+'-'+mjd+'.fits')
    if not ptt.exists(spZline_file):
        if ptt.exists(spZline_file+'.gz'):
            spZline_file = spZline_file+'.gz'
        else:
            splog.log('ERROR: Missing '+spZline_file)
            exit()
    spZline = Table(fits.getdata(spZline_file,1))
    spZline_hdr = fits.getheader(spZline_file,1)
    
    return(spZline, spZline_hdr)

def read_spZbest(sp1d_dir, field, mjd):
    spZbest_file = ptt.join(sp1d_dir, 'spZbest-'+field+'-'+mjd+'.fits')
    if not ptt.exists(spZbest_file):
        if ptt.exists(spZbest_file+'.gz'):
            spZbest_file = spZbest_file+'.gz'
        else:
            splog.log('ERROR: Missing '+spZbest_file)
            exit()
    spZbest = fits.getdata(spZbest_file,2)
    spZbest_hdr = fits.getheader(spZbest_file,2)
    spZbest_map = fits.getdata(spZbest_file,1)

    return(spZbest, spZbest_map, spZbest_hdr)

def set_description(card, h, field_hdr):

    try:
        if field_hdr[h] == card:
            return(field_hdr.cards[h])
    except:
        pass
        
    added_card = {'TILE':"Tile ID for SDSS BOSS plates (-1 for SDSS)",
                  'THETA_COVAR':"Covariance matrix for THETA",
                  'FLUX':"coadded calibrated flux",
                  'LOGLAM':"log10(wavelength [Angstrom])",
                  'IVAR':"inverse variance of flux",
                  'AND_MASK':"AND mask",
                  'OR_MASK':"OR mask",
                  'WDISP':"Wavelength dispersion in number of pixel",
                  'SKY':"subtracted sky flux",
                  'MODEL':"pipeline best model fit for class and z",
                  'WRESL':"spectral resolution in A units",
                  }

    try:
        match = np.where(str(card) == np.array(field_hdr.cards)[:,1])[0]
    except:
        match = []
    if len(match) > 0:
        col = field_hdr.cards[match[0]][1]
        comment = field_hdr.cards[match[0]][2]
        return(h,col,comment)
    elif card in added_card.keys():
        return(h,card,added_card[card])
    else:
        return(h, card,None)
        
def spSpec_reformat(boss_spectro_redux, run2d, run1d, field, mjd,
                    plot=False, epoch=False, lsdr10=False,
                    allsky=False, custom=None):
    if allsky is False:
        field = field_to_string(field)

    specfull_dir = field_spec_dir(boss_spectro_redux, run2d, field, mjd,
                                  epoch=epoch, custom=allsky, custom_name=custom)
    speclite_dir = field_spec_dir(boss_spectro_redux, run2d, field, mjd,
                                  epoch=epoch, full=False, custom=allsky,
                                  custom_name=custom)
    specImg_dir  = field_png_dir(boss_spectro_redux, run2d, run1d, field, mjd,
                                 epoch=epoch, custom=allsky, custom_name=custom)
    fdir = field_dir(ptt.join(boss_spectro_redux, run2d), field, custom=allsky)
    if epoch is True:
        sp1d_dir =  ptt.join(fdir, 'epoch', run1d)
        logfile = ptt.join(fdir, 'epoch', 'spSpec_reformat-'+field+'-'+mjd+'.log')
        sfiles = glob(ptt.join(fdir,'epoch','coadd',mjd,'spSpec-'+field+'-'+mjd+'-*.fits'))
    else:
        sp1d_dir =  ptt.join(fdir, run1d)
        logfile = ptt.join(fdir, 'spSpec_reformat-'+field+'-'+mjd+'.log')
        sfiles = glob(ptt.join(fdir,'coadd',mjd,'spSpec-'+field+'-'+mjd+'-*.fits'))

    splog.open(logfile = logfile)
    splog.log('Log file '+logfile+' opened '+ time.ctime())
    
    spAll, spAll_hdr, SDSSC2BV = read_spAll(specfull_dir, field, mjd)
    spZall, spZall_hdr = read_spZall(sp1d_dir, field, mjd)
    spZline, spZline_hdr = read_spZline(specfull_dir, field, mjd)
    spZbest, spZbest_map, spZbest_hdr = read_spZbest(sp1d_dir, field, mjd)

    model_targIdx = np.asarray(spZbest_map['TARGET_INDEX'].data)
    spAll_targ_hdr = spAll_hdr
    spZall_targ_hdr = spZall_hdr
    spZline_targ_hdr = spZline_hdr
    
    makedirs(specfull_dir, exist_ok=True)
    makedirs(speclite_dir, exist_ok=True)
    makedirs(specImg_dir,  exist_ok=True)

    units = {'FLUX':"10^-17 ergs/s/cm^2/Angs",
             'LOGLAM':"log10(Angs)",
             'WDISP':"Pixels",
             'SKY':"10^-17 ergs/s/cm^2/Angs",
             'WRESL':"Angs"}
             
    if len(sfiles) == 0:
        splog.log('ERROR: No files matching spSpec-'+field+'-'+mjd+'-*.fits')
        exit()
    files = []
    
    files = Table(names=('name', 'RA', 'DEC'), dtype=(str, float, float))

    for i, spSpecF in enumerate(sfiles):
        specF = ptt.basename(spSpecF).replace('spSpec','spec')
        splog.info(f'{specF} ({i+1}/{len(sfiles)})')
        with fits.open(spSpecF) as spSpec:
            spSpec.verify('silentfix')
            spec = [spSpec[0]]#,spSpec[1]]
            spec[0].header.set('SDSSC2BV', SDSSC2BV, 'SDSS5_TARGET_FLAG Carton to Bit Version')
            
            #Coadd Extension
            cols = []
            cols.append(fits.Column(name='FLUX',     format='E', array=spSpec[1].data['FLUX'],      unit="10^-17 ergs/s/cm^2/Angs") )
            cols.append(fits.Column(name='LOGLAM',   format='E', array=spSpec[1].data['LOGLAM'],    unit="log10(Angs)") )
            cols.append(fits.Column(name='IVAR',     format='E', array=spSpec[1].data['IVAR'] ))
            cols.append(fits.Column(name='AND_MASK', format='J', array=spSpec[1].data['AND_MASK']) )
            cols.append(fits.Column(name='OR_MASK',  format='J', array=spSpec[1].data['OR_MASK']) )
            cols.append(fits.Column(name='WDISP',    format='E', array=spSpec[1].data['WDISP'],     unit="Pixels") )
            cols.append(fits.Column(name='SKY',      format='E', array=spSpec[1].data['SKY'],       unit="10^-17 ergs/s/cm^2/Angs") )
            idx = np.where(model_targIdx == spSpec[2].data['TARGET_INDEX'][0])[0]
            if len(idx) == 0:
                model = np.zeros(len(spSpec[1].data['FLUX']),dtype=float)
            else:
                model = spZbest[idx[0]]
            cols.append(fits.Column(name='MODEL',    format='E', array=model) )
            cols.append(fits.Column(name='WRESL',    format='E', array=spSpec[1].data['WRESL'],     unit="Angs") )
            spec.append(fits.BinTableHDU.from_columns(cols, name = 'COADD'))
            for h in spec[-1].header:
                if 'TTYPE' in h:
                    card = set_description(spec[-1].header.get(h), h, None)
                    spec[-1].header.set(card[0],card[1],card[2])

            #SPALL Extension
            spAll_target = spAll[spAll['TARGET_INDEX'].data == spSpec[2].data['TARGET_INDEX'][0]]
            spAll_targ_hdr['NAXIS2'] =1
            spec.append(fits.BinTableHDU(spAll_target, header=spAll_targ_hdr, name = 'SPALL'))
            for h in spAll_targ_hdr:
                if 'TTYPE' in h:
                    card = set_description(spAll_targ_hdr.get(h), h, spAll_targ_hdr)
                    spec[-1].header.set(card[0],card[1],card[2])
            spec[-1].header.set('SDSSC2BV', SDSSC2BV, 'SDSS5_TARGET_FLAG Carton to Bit Version', before='EXTNAME')
            spec[-1].add_checksum()
            
            #ZALL Extension
            spZall_target = spZall[spZall['TARGET_INDEX'] == spSpec[2].data['TARGET_INDEX'][0]]
            spZall_targ_hdr['NAXIS2'] = len(spZall_target)
            spec.append(fits.BinTableHDU(spZall_target, header=spZall_targ_hdr, name = 'ZALL'))
            for h in spZall_targ_hdr:
                if 'TTYPE' in h:
                    card = set_description(spZall_targ_hdr.get(h), h, spAll_targ_hdr)
                    if card[1] is None: continue
                    spec[-1].header.set(card[0],card[1],card[2])
            spec[-1].header.remove('COMMENT', ignore_missing=True, remove_all=True)
            spec[-1].add_checksum()

            #ZLINE Extension
            spZline_target = spZline[spZline['TARGET_INDEX'] == spSpec[2].data['TARGET_INDEX'][0]]
            spZline_targ_hdr['NAXIS2'] = len(spZline_target)
            spec.append(fits.BinTableHDU(spZline_target, header=spZline_targ_hdr, name = 'ZLINE'))
            for h in spZline_targ_hdr:
                if 'TTYPE' in h:
                    card = set_description(spZline_targ_hdr.get(h), h, spZline_targ_hdr)
                    if card[1] is None: continue
                    spec[-1].header.set(card[0],card[1],card[2])
            spec[-1].add_checksum()
            
            last_ext  = len(spec)
            # Exposure Extensions
            hdul = fits.HDUList(spec)
            hdul.verify('silentfix')
            hdul.writeto(ptt.join(speclite_dir,specF), overwrite=True, checksum=True)
            spec.extend(spSpec[3:])
            

            for i, hdu in enumerate(spec[last_ext:]):
                for h in hdu.header:
                    if 'TTYPE' in h:
                        card = set_description(hdu.header.get(h), h, None)
                        hdu.header.set(card[0],card[1],card[2])
                        if card[1] in units:
                            hdu.header.set(card[0].replace('TYPE','UNIT'),units[card[1]],card[1])

            for i, hdu in enumerate(spec):
                spec[i].add_checksum()
            hdul = fits.HDUList(spec)
            hdul.verify('silentfix')
            hdul.writeto(ptt.join(specfull_dir,specF), overwrite=True, checksum=True)
            
            catid = '-'.join(specF.replace('.fits','').replace('.gz','').split('-')[3:])
            
            if plot:
                try:
                    files = SDSS_specplot(specImg_dir, Table(spec[1].data), Table(spec[2].data)[0],catid,
                                      hdr = spec[0].header, files = files, allsky = allsky, field=field)
                except:
                    print(specF, spSpec[2].data['TARGET_INDEX'][0])
                    print(spAll['TARGET_INDEX'].data)
                    #exit()
                    
    if plot:
        build_html(specImg_dir, field, mjd, files=files, lsdr10=lsdr10, allsky=allsky)
    
    splog.log('Successful completion of spSpec_reformat at '+ time.ctime())
    splog.close()
     
def build_html(specImg_dir, field, mjd, files=Table(), lsdr10=False, allsky=False):
    if allsky is False:
        pmjd = field_to_string(field)+'-'+mjd
    else:
        pmjd = field+'-'+mjd
    with open(ptt.join(specImg_dir, 'tmp-index.html'), 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>'+'\n')
        f.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'+'\n')
        f.write('<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">'+'\n')
        f.write('<head>'+'\n')
        f.write('<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />'+'\n')
        f.write('<title>'+pmjd+'</title>'+'\n')
        f.write('<style type="text/css">'+'\n')
        f.write('body { background: #111; }'+'\n')
        f.write('td { color: gray; }'+'\n')
        f.write('</style>'+'\n')
        f.write('</head>'+'\n')
        f.write('<body>'+'\n')
        f.write('<table border="0" cellspacing="5">'+'\n')

        nper=5
        for i, row in enumerate(files):
            if((i % nper) == 0): f.write('<tr>'+'\n')
                        
            stamp = '.'
            
            stampLS='https://www.legacysurvey.org/viewer?ra='+'%.5f' % row['RA']+'&dec='+ '%.5f' % row['DEC']
            stampLS+='&layer=ls-dr10&spectra&mark='+'%.5f' % row['RA']+','+'%.5f' % row['DEC']+'&zoom=14'
            
            stamp='http://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg?TaskName=Skyserver.Chart.Image&ra='
            stamp+='%.5f' % row['RA'] +'&dec='+ '%.5f' % row['DEC']
            stamp+='&scale=0.1&width=512&height=512&opt=G&query=&Grid=on'
            currbase=row['name']
            outbase=specImg_dir+'/'+currbase
            pmjdf = row['name'].replace('spec-image-','')
            if lsdr10:
                f.write('<td><a href="'+stamp+'">'+pmjdf+ '</a> <br/><a href="'+stampLS+'">(LS-DR10)</a><br /><a href="'+currbase+'.png">'+'\n')
            else:
                f.write('<td><a href="'+stamp+'">'+pmjdf+ '</a> <br /><a href="'+currbase+'.png">'+'\n')
            f.write('<img src="'+currbase+'.thumb.png" alt="'+pmjdf+'" /></a>'+'\n')
            f.write('</td>'+'\n')
            if(((i % nper) == nper-1) or (i == len(files)-1)): f.write('</tr>'+'\n')
        f.write('</table>'+'\n')
        f.write('</body>'+'\n')
        f.write('</html>'+'\n')
    rename(ptt.join(specImg_dir, 'tmp-index.html'), ptt.join(specImg_dir, 'index.html'))


def build_title(spAll, catid, allsky=False, field=None):

    survey = spAll['SURVEY'].strip().replace('_','\_')
    program = spAll['PROGRAMNAME'].strip().replace('_','\_')

    if 'SDSS_ID' in spAll.colnames:
        SDSSID = spAll['SDSS_ID']
        if ((SDSSID == -999) or (SDSSID == -1) or (SDSSID == 0)): SDSSID = ''
    else:
        SDSSID = None

    if allsky is False:
        field = spAll['FIELD']
        legacy = False if int(field)>=15000 else True
        plates = False if int(field)>=16000 else True
        sfield = field_to_string(spAll['FIELD'])
    else:
        legacy = False
        plates = False
        field  = field
        sfield = field
    #ptitle = r'Survey: $'+survey+r'$, Program: $'+program+r'$, '
    
    try:
        obs =spAll['OBS'].strip()
    except:
        obs = ''
    obstitle = r'Observatory: '+obs+r', '
    ptitle = r'Program: $'+program+r'$, '

    if survey == 'sdss':
        targ_title='Target'
        primtarget = spAll['PRIMTARGET']
        targets= ' '.join(sdss_flagname('TARGET', primtarget))
    elif survey == 'segue1':
        targ_title='Target'
        primtarget = spAll['PRIMTARGET']
        targets= ' '.join(sdss_flagname('SEGUE1_TARGET', primtarget))
    elif survey ==  'boss':
        targ_title='Target'
        boss_target1 = spAll['BOSS_TARGET1']
        ANCILLARY_TARGET1 = spAll['ANCILLARY_TARGET1']
        targets= ' '.join(sdss_flagname('BOSS_TARGET1', boss_target1))+ ' '+ ' '.join(sdss_flagname('ANCILLARY_TARGET1', ANCILLARY_TARGET1))
    elif ((survey.lower() in ['bhm-mwm', 'bhm', 'mwm', 'mwm-bhm','open\_fiber'])):
        targ_title='Firstcarton'
        targets= spAll['FIRSTCARTON'].replace('_','\_')
    else:
        targ_title = ''
        targets = 'NA'
    if targets.strip() == 'NA':
        targ_title='Firstcarton'
        targets=str(spAll['OBJTYPE']).strip().replace('_','\_')

    ptitle= ptitle+targ_title+r': $'+targets+r'$'

    if not legacy:
        if SDSSID is not None:
            ptitle = ptitle+ '\n '
            ptitle = ptitle+ 'SDSS_ID='+str(SDSSID)
        ptitle = ptitle+', '+'CatID='+str(catid)+ '\n '
    else:
        if SDSSID is not None:
            ptitle = ptitle+', '+'SDSS_ID='+str(SDSSID)+ '\n '
        else:
            ptitle = ptitle+ '\n '

    if 'FIBER_RA' in spAll.columns:
        title1= ('RA='+"%.5f" % spAll['FIBER_RA']+', '+'Dec='+"%.5f" % spAll['FIBER_DEC']+', '+
                 'Field='+sfield+', '+obstitle+'TargetIndex='+str(spAll['TARGET_INDEX'])+', '+
                 'MJD='+str(spAll['MJD']))
    else:
        title1= ('RA='+"%.5f" % spAll['PLUG_RA']+', '+'Dec='+"%.5f" % spAll['PLUG_DEC']+', '+
                 'Field='+sfield+', '+obstitle+'Fiber='+str(spAll['FIBERID'])+', '+
                 'MJD='+str(spAll['MJD']))
    ptitle = ptitle + title1+'\n'

    if(spAll['Z'] < 1000./299792.):
        zstr= r'$cz=$'+str(int(spAll['Z']*299792.))+r'$\pm$'+str(int(spAll['Z_ERR']*299792.))+' km/s'
    else:
        zstr= r'$z=$'+"%.5f" % spAll['Z']+r'$\pm$'+"%.5f" % spAll['Z_ERR']
        
    title2 = zstr+', Class='+spAll['CLASS'].strip()+' '+spAll['SUBCLASS'].strip()
    if len(spAll['SUBCLASS'].strip()) == 0:
        title2 = zstr+', Class='+spAll['CLASS'].strip()

    
    mag_vec=spAll['FIBER2MAG']
    m_i=mag_vec[3]
    title2=title2+r', mag$_{i,fib2}$='+"%.2f" % m_i
    
    sn = spAll['SN_MEDIAN_ALL']
    title2=title2+r', S/N='+"%.2f"%sn
    
    ptitle = ptitle + title2+'\n'

    if spAll['ZWARNING'] > 0:
        warnings= ' '.join(sdss_flagname('ZWARNING', spAll['ZWARNING']))
        title3= 'Warnings: '+warnings
    else:
        title3='No warnings'
    title3 = title3+', idlspec2d='+spAll['RUN2D']
    ptitle = ptitle + title3
    
    ptitle = ptitle.replace(r'$$',r'')
    return(ptitle)
    

def SDSS_specplot(basedir, Coadd_Table, spAll, catalogID, files = Table(), xra=[3501., 10499.],
                  hdr=None, allsky=False, field=None):
    if not allsky:
        outbase = 'spec-image-'+field_to_string(spAll['FIELD'])+'-'+str(spAll['MJD'])+'-'+catalogID
    else:
        outbase = 'spec-image-'+field+'-'+str(spAll['MJD'])+'-'+catalogID

    flux = Coadd_Table['FLUX'].data
    loglam = Coadd_Table['LOGLAM'].data
    ivar = Coadd_Table['IVAR'].data
    wave = np.power(10.0, loglam)
    igd = np.where(ivar > 0)[0]
    if len(igd)> 0:
        ist = min(igd)
        ind = max(igd)
        flux = flux[ist:ind]
        ivar = ivar[ist:ind]
        wave = wave[ist:ind]
    sflux = sdss_spec_smooth(np.log10(wave), flux, 100)
    
    sscale=1.5
    fig, axs = plt.subplots(1, figsize=(10.5*sscale, 7.5*sscale), dpi=72*2)
    
    axs.plot(wave, sflux, color='k', alpha=1, zorder = 2, lw=1)
    
    igd= np.where((ivar > 0) & (np.abs(wave-5577.) > 4.) & (wave < 10000.) & (wave > 3700.))[0]
    if len(igd) > 0:
        yra= [np.nanmin(sflux[igd]),np.nanmax(sflux[igd])]
    else:
        yra= [np.nanmin(sflux),np.nanmax(sflux)]
    size= 0.07*(yra[1]-yra[0])
    yra=yra+np.array([-1.2,1.7])*size*1.7
    if(yra[0] < -2.):
        yra[0]=-1.999
    if yra[0] == yra[1]:
        yra=[0,1]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        err= np.sqrt(1/ivar)
        ibad= np.where(ivar <= 0)[0]
        nans, x= nan_helper(err)

        if (len(nans) > 0) and (sum(nans) != len(nans)):
            err[nans]= np.interp(x(nans), x(~nans), err[~nans])

        err=sdss_spec_smooth(np.log10(wave), err, 400)
        yerr_u = err.copy()
        yerr_l = err.copy()
        yerr_u[ibad] = yra[1] - sflux[ibad]
        yerr_l[ibad] = sflux[ibad] - yra[0]

    axs.fill_between(rebin(wave,2), rebin(sflux-yerr_l,2), rebin(sflux+yerr_u,2), color='k', alpha=.2, zorder = 1, ec=None, linewidth=0., step='mid')
    axs.fill_between(wave, sflux-yerr_l,sflux+yerr_u, color='k', alpha=.2, zorder = 1, ec=None, linewidth=0.)
    alines = {4300.:"G",5895.:"Na D",5175.:"Mg",8498.:"CaII",8542.:"",8662.:"",3968.:" H",3938.:"K "}
    
    elines= {3727.:'OII', 3869.7867:'NeIII', 4105.8884: r'H$\mathbf{\delta}$',
             4341.6803: r'H$\mathbf{\gamma}$', 4364.3782:'OIII', 4862.6778: r'H$\mathbf{\beta}$',
             4960.2140:'', 5008.1666:'OIII', 5876.:'HeI', 6301.9425:'OI',
             6549.7689:'NII', 6564.6127:r'H$\mathbf{\alpha}$', 6585.1583:'NII',
             6718.1642: '', 6732.5382: 'SII',7137.6370:'ArIII', 2800.:'MgII',
             1216.:r'Ly$\mathbf{\alpha}$', 1549.:'CIV', 1640.:'HeII', 1909.:'CIII',
             2326.:'CII', 1400.:'SiIV+OIV'}

    ewaves = np.asarray(list(elines.keys()))
    elines = np.asarray(list(elines.values()))
    xsize= 0.07*(xra[1]-xra[0])
    
    yoff= np.zeros_like(ewaves)
    xoff= np.zeros_like(ewaves)

    ioff= np.where(np.abs(ewaves-6585.) < 3.)[0]
    xoff[ioff]=0.3*xsize
    yoff[ioff]=-0.5*size
    ioff= np.where(np.abs(ewaves-6549.) < 3.)[0]
    xoff[ioff]=-0.3*xsize
    yoff[ioff]=-0.5*size
    ioff= np.where(np.abs(ewaves-4862.) < 3.)[0]
    xoff[ioff]=-0.1*xsize
    ioff= np.where(np.abs(ewaves-5008.) < 3.)[0]
    xoff[ioff]=0.05*xsize
    ioff= np.where(np.abs(ewaves-4364.) < 3.)[0]
    xoff[ioff]=0.20*xsize
    ioff= np.where(np.abs(ewaves-4341.) < 3.)[0]
    xoff[ioff]=-0.05*xsize
    yoff[ioff]=size

    ewave=ewaves*(1.+spAll['Z'])
    iwave= np.where((ewave > xra[0]) & (ewave < xra[1]))[0]
    for i in iwave:
        inear= np.where((wave > ewave[i]-100.) & (wave < ewave[i]+100) &
                        (ivar > 0) & (np.abs(wave-5577.) > 4.))[0]
        if len(inear) > 0 :
            e_val= np.nanmax(sflux[inear])
            axs.annotate(elines[i],
                        xy= (ewave[i], e_val*1.05+yoff[i]), xycoords='data',
                        xytext=(ewave[i], e_val*1.05+yoff[i]+size), textcoords='data',
                        arrowprops={'arrowstyle':'-', 'color':'b', 'linewidth':1.5},
                        color='k',horizontalalignment='center', fontfamily='serif',
                        clip_on=True,annotation_clip=True,fontweight= 'heavy',
                        )

    awaves = np.asarray(list(alines.keys()))
    alines = np.asarray(list(alines.values()))
    
    yoff = np.full_like(awaves, -1*size*0.2)
    xoff = np.zeros_like(awaves)
    ioff= np.where(np.abs(awaves-3938.) < 3.)[0]
    xoff[ioff]=-0.07*xsize
    ioff= np.where(np.abs(awaves-3968.) < 3.)[0]
    xoff[ioff]=0.07*xsize

    awave=awaves*(1.+spAll['Z'])
    iwave= np.where((awave > xra[0]+100.) & (awave < xra[1]-100.))[0]
    for i in iwave:
        inear= np.where((wave > awave[i]-100.) & (wave < awave[i]+100) &
                        (ivar > 0) & (np.abs(wave-5577.) > 4.))[0]
        if len(inear) > 0:
            aval= np.nanmin(sflux[inear])
            ylow = aval*0.95+yoff[i]-size
            ann_pars = {'xy':(awave[i], aval*0.95+yoff[i]), 'xycoords':'data',
                        'xytext':(awave[i], ylow), 'textcoords':'data',
                        'arrowprops':{'arrowstyle':'-', 'color':'r', 'linewidth':1.5},
                        'color':'k','horizontalalignment':'center', 'fontfamily':'serif',
                        'clip_on':True,'annotation_clip':True, 'fontweight':'heavy'}
            if ylow<yra[0]:
                ylow = yra[0]
                alines[i] = ''
                ann_pars['fontsize']= 0
                ann_pars['xy']=(awave[i], aval*0.95+yoff[i])
                ann_pars['xytext']=(awave[i], ylow)
            axs.annotate(alines[i], **ann_pars)

    axs.plot(wave, sflux, color='k', alpha=1, zorder = 2, lw=1)


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        axs.set_ylim(yra)
    axs.set_xlim(xra)

    axs.set_title(build_title(spAll,catalogID, allsky=allsky, field=field), horizontalalignment='center')
    axs.set_ylabel(r'$f_\lambda~(10^{-17}~ergs/s/cm^2/\mathrm{\AA})$')
    axs.set_xlabel(r'Wavelength $(\mathrm{\AA})$')
    fig.tight_layout()
    hmargin = max([fig.subplotpars.left, 1.0-fig.subplotpars.right])
    vmargin = max([fig.subplotpars.bottom, 1.0-fig.subplotpars.top])
    plt.subplots_adjust(left = hmargin, right = (1.0-hmargin), top= (1.0-vmargin), bottom = vmargin)
    plt.savefig(ptt.join(basedir,outbase+'.png'))
    plt.close(fig)

    
    
    fig = image.thumbnail(ptt.join(basedir,outbase+'.png'), ptt.join(basedir, outbase+'.thumb.png'), scale=0.08)
    plt.close(fig)


    if hdr is not None:
        info = PngImagePlugin.PngInfo()
        for key in hdr.keys():
            if key in ['CHECKSUM','DATASUM','SIMPLE','BITPIX','NAXIS','EXTEND','']:
                continue
            info.add_text(key,str(hdr[key]))
        im = Image.open(ptt.join(basedir,outbase+'.png'))
        im.save(ptt.join(basedir,outbase+'.png'), "PNG", pnginfo=info)
 
 
    if 'FIBER_RA' in spAll.columns:
        files.add_row((outbase, spAll['FIBER_RA'], spAll['FIBER_DEC']))
    else:
        files.add_row((outbase, spAll['PLUG_RA'], spAll['PLUG_DEC']))

    
    return(files)

