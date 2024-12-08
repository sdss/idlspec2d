#
from boss_drp.utils import find_nearest_indx
from boss_drp.field import field_to_string

from pydl.pydlutils.trace import traceset2xy, TraceSet

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os.path as ptt
import os
import glob
import warnings
import numpy as np

def bc_smooth(data, window_size):
    if window_size % 2 == 0:
        window_size += 1
    
    # Create a boxcar (flat) window
    window = np.ones(window_size) / window_size
    
    # Apply the convolution
    smoothed_data = np.convolve(data, window, mode='same')
    
    return smoothed_data


def plot(mjd, expid, ccd, redo=False, outdir = '/data/boss/sos/tests/', mask_end=False):
    mjd = str(mjd)
    sos_dir = os.getenv('BOSS_SOS_S') if '2' in ccd else os.getenv('BOSS_SOS_N')
    if sos_dir is None:
        sos_dir = '/data/boss/sos'
        if redo:
            sos_dir = '/data/boss/sosredo'
    sci_f = ptt.join(sos_dir,mjd,f'sci-*-{ccd}-{str(expid).zfill(8)}.fits')
    sci_f = glob.glob(sci_f)[0]
    hdr = fits.getheader(sci_f)
    confid = hdr['CONFID']
    fieldid = field_to_string(hdr['FIELDID'])
    MJD = hdr['MJD']

    arc_f = glob.glob(ptt.join(sos_dir, mjd,f'wset-{MJD}-{fieldid}-*-{ccd}.fits'))
    if len(arc_f) == 0:
        arc_f = glob.glob(ptt.join(sos_dir, mjd,f'wset-{MJD}-*-*-{ccd}.fits'))
        exps = np.array([ptt.basename(x).split('-')[0] for x in arc_f]).astype(int)
        idx = find_nearest_indx(exps, int(expid))
        arc_f = arc_f[idx]
    else:
        arc_f = arc_f[0]
    wset = fits.getdata(arc_f,1)
    xx, loglam = traceset2xy(TraceSet(wset))
    
    fmap_f = ptt.join(sos_dir, mjd, f'spfibermap-{fieldid}-{mjd}-{ccd}.fits')
    try:
        fmap = fits.getdata(fmap_f,f'CONFSUMMARYF-{confid}.PAR')
    except:
        try:
            fmap = fits.getdata(fmap_f,f'CONFSUMMARY-{confid}.PAR')
        except:
            fmap = fits.getdata(fmap_f,2)

    
    sci = fits.getdata(sci_f,0)
    sn = fits.getdata(sci_f,2)
    ivar = fits.getdata(sci_f,1)
    
    num_panels = 500
    panels_per_page = 6
    num_pages = (num_panels + panels_per_page - 1) // panels_per_page  # Calculate total number of pages

    # Create a PDF file
    os.makedirs(outdir,exist_ok=True)
    with PdfPages(ptt.join(outdir,f'{expid}-{ccd}-{mjd}.pdf')) as pdf:
        fiber = 0
        for page in range(num_pages):
            fig, axes = plt.subplots(6, 1, figsize=(8.5,11))  # 2x4 grid, each subplot on A4 paper
            axes = axes.flatten()  # Flatten to 1D array for easier iteration

            for i in range(panels_per_page):
                panel_idx = page * panels_per_page + i
                if panel_idx < num_panels:
                    ax = axes[i]
                    data = sci[fiber]
                    sm_data = bc_smooth(data, 10)
                    wave = np.power(10,loglam[fiber])
                    mask = np.where((wave < 10000) & (wave> 5000))[0]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        ax.plot(wave, np.sqrt(1./ivar[fiber]), color='r',alpha=.2)
                    ax.plot(wave, data, color='k', alpha=.2)
                    ax.plot(wave, sm_data, color='k')

                    title = f"FiberID={fiber+1} catid={fmap[fiber]['CATALOGID']} Carton={fmap[fiber]['FIRSTCARTON']}\nmag={','.join(np.asarray(fmap[fiber]['CATDB_MAG']).astype(str).tolist())} sn={sn[fiber]}"
                    
                    ax.set_title(title)
                    if mask_end:
                        sm_data = sm_data[mask]
                    ax.set_ylim(-1,np.nanmean(sm_data[sm_data >= 0])+1.5*np.nanstd(sm_data[sm_data >= 0]))
                    fiber += 1
                else:
                    axes[i].axis('off')  # Hide unused subplots
            plt.tight_layout(pad=1.0)
            pdf.savefig(fig)
            plt.close(fig)
    print('Saved to '+ptt.join(outdir,f'{expid}-{ccd}-{mjd}.pdf'))

