#
from boss_drp.utils.lock import lock, unlock
from boss_drp.field import field_to_string
from boss_drp.utils.splog import splog
import boss_drp
from boss_drp import idlspec2d_dir, favicon


from pydl.pydlutils.trace import traceset2xy, TraceSet

from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['pdf.fonttype'] = 42

from matplotlib.backends.backend_pdf import PdfPages
import os.path as ptt
import os
import glob
import warnings
import numpy as np
from datetime import datetime
import warnings
from jinja2 import Template


class SciFrame:
    # This class encapsulates the data and metadata for a single SOS science frame, 
    # including the wavelength solution, fiber mapping, and smoothed spectra.
    def __init__(self, mjd, ccd, sci_f, sos_dir, window_size=10): 
        # Read the science frame header to extract metadata and determine file paths for the arc and fiber mapping files
        # Also read the science data, inverse variance, and signal-to-noise arrays, and set up a smoothing window for later use in plotting
        # 
        # Parameters:
        # -----------
        # mjd: str
        #     The MJD of the reduction, used to locate the relevant files in the SOS directory structure.
        # ccd: str
        #     The CCD identifier ('b' or 'r') to specify which science frame to read.
        # sci_f: str
        #     The file path to the science frame FITS file to read.
        # sos_dir: str
        #     The base directory for SOS reductions, used to locate the arc and fiber mapping files.
        # window_size: int, optional
        #     The size of the smoothing window to apply to the spectra for plotting (default is 10).                

        hdr = fits.getheader(sci_f)
        self.sci_f = sci_f
        self.ccd = ccd
        self.confid = hdr['CONFID']
        self.fieldid = field_to_string(hdr['FIELDID'])

        arc_f = ptt.join(sos_dir,hdr['WSETFILE'])
        wset = fits.getdata(arc_f,1)
        _, self.loglam = traceset2xy(TraceSet(wset))
        self.wave = np.power(10,self.loglam)
        
        fmap_f = ptt.join(sos_dir, mjd, f'spfibermap-{self.fieldid}-{mjd}-??.fits')
        fmap_f = glob.glob(fmap_f)[0]
        try:
            self.fmap = fits.getdata(fmap_f,f'CONFSUMMARYF-{self.confid}.PAR')
        except:
            try:
                self.fmap = fits.getdata(fmap_f,f'CONFSUMMARY-{self.confid}.PAR')
            except:
                self.fmap = fits.getdata(fmap_f,2)
        self.fmap = self.fmap[self.fmap['FIBERTYPE'] == 'BOSS']
        
        self.sci = fits.getdata(sci_f,4)
        self.sn = fits.getdata(sci_f,2)
        self.ivar = fits.getdata(sci_f,5)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.err = np.sqrt(1./self.ivar)
        self.run2d = hdr.get('RUN2D',boss_drp.__version__)

        if window_size % 2 == 0:
            window_size += 1
        self._smooth_window = np.ones(window_size) / window_size

    def get_fiber(self, fiber):
        # Given a fiber index, this method retrieves the corresponding spectrum, wavelength solution,
        # error array, fiber mapping information, and signal-to-noise ratio for that fiber from the science frame data. 
        # It also applies a smoothing operation to the spectrum using a predefined window size for later use in plotting.
        #
        # Parameters:
        # -----------
        # fiber: int
        #     The index of the fiber to retrieve data for
        #
        data = self.sci[fiber]
        sm_data = np.convolve(data, self._smooth_window, mode='same')
        wave = self.wave[fiber]
        return wave, data, sm_data, self.err[fiber], self.fmap[fiber], self.sn[fiber]

colors = {'b':"#0072B2",'r':"#CC79A7",'b_err':"#E69F00",'r_err':"#F0E442"}

def _plotone(ax, spectra_dict, fiberids, fiber, mask_end=False):
    # This function takes a matplotlib axis, a dictionary of SciFrame objects for each CCD,
    # a list of fiber IDs, and an index for the current fiber to plot.
    # It retrieves the spectrum, wavelength solution, error array, fiber mapping information,
    # and signal-to-noise ratio for the specified fiber from each CCD's SciFrame, applies any necessary scaling to
    # ensure consistent flux levels across CCDs, and plots the raw and smoothed spectra along with the error array on the provided axis.
    # The function also sets the plot title to include relevant information about the fiber and its associated metadata, 
    # and adjusts the plot limits based on the data being plotted.
    # 
    # Parameters:
    # -----------
    # ax: matplotlib.axes.Axes
    #     The matplotlib axis on which to plot the spectra for the specified fiber.
    # spectra_dict: dict
    #     A dictionary where the keys are CCD identifiers ('b' or 'r') and the values are SciFrame objects containing the data 
    #       and metadata for each CCD's science frame.
    # fiberids: list or array-like
    #     A list or array of fiber IDs that are being plotted, used to index into the SciFrame data for each CCD.
    # fiber: int
    #     The index of the current fiber within the fiberids list to plot.
    # mask_end: bool, optional
    #     If True, the function will mask out the ends of the spectra (outside of the wavelength range of 5000-10000 Angstroms)
    #        during plotting to focus on the central portion of the spectra (default is False).
    ymin = []
    ymax = []
    scale = None
    for ccd in spectra_dict:
        wave, data, sm_data, err, fmap, sn = spectra_dict[ccd].get_fiber(fiberids[fiber])
        mask = (sm_data >= 0) & (wave > 6000) & (wave < 6500)
        if np.any(mask):
            if scale is None:
                scale = np.nanmean(sm_data[mask]) 
            else:
                scale = scale/np.nanmean(sm_data[mask]) 
                data = data * scale
                sm_data = sm_data * scale
                err = err * scale
        else:
            scale = None


        mask = np.where((wave < 10000) & (wave> 5000))[0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ax.plot(wave, err, color=colors[f'{ccd}_err'],alpha=.2)
        if mask_end:
            data = data[mask]
            sm_data = sm_data[mask]
            wave = wave[mask]
        ax.plot(wave, data, color=colors[ccd], alpha=.2)
        ax.plot(wave, sm_data, color=colors[ccd])
        ymin.append(-1)
        if ccd == 'b':
            sm_data = sm_data[(wave > 3600) & (wave < 6500)]
        else:
            sm_data = sm_data[(wave > 6000) & (wave < 10000)]
        
        ymax.append(np.nanmean(sm_data[sm_data >= 0])+1.5*np.nanstd(sm_data[sm_data >= 0]))


    mags = np.asarray(fmap['CATDB_MAG'], dtype=float)
    mags = ','.join('nan' if np.isnan(v) else f"{v:.2f}" for v in mags)
    title = f"FiberID={fiberids[fiber]+1} catid={fmap['CATALOGID']} Carton={fmap['FIRSTCARTON']}\nmag={mags} sn={sn}"
    ax.set_title(title)
    if mask_end:
        sm_data = sm_data[mask]
    ax.set_xlim(3600,10000)
    ax.set_ylim(min(ymin),max(ymax))

    return (fmap['RA'], fmap['DEC'], fmap['CATALOGID'])


def plot(mjd, expid, obs, ccd_name, sos_dir= None, redo=False, mask_end=False, ToOs= False, assigned=False, science=False, pdf=True, single_ccd=False):
    # This function serves as the main entry point for plotting the spectra for a given MJD, exposure ID, observatory, and CCD.
    # It locates the relevant science frame files for the specified parameters, creates SciFrame objects
    # to encapsulate the data and metadata for each CCD, applies any specified filtering criteria to select which fibers to plot,
    # and then generates either individual PNG plots for each fiber or a multi-panel PDF containing all
    # the selected fibers, depending on the user's choice. The function also handles file locking to ensure that plots are saved safely without conflicts.
    # Parameters:
    # -----------
    # mjd: str
    #     The MJD of the reduction, used to locate the relevant files in the SOS    directory structure.
    # expid: str or int
    #     The exposure ID to plot, used to identify the specific science frame files to read and plot.
    # obs: str
    #     The observatory ('LCO' or 'APO') to determine which CCDs to plot and where to find the relevant files.
    # ccd_name: str
    #     The CCD identifier ('b' or 'r') to specify which science frame to read and plot.
    # sos_dir: str, optional
    #     The base directory for SOS reductions, used to locate the relevant files 
    #       (default is None, in which case it will be determined from environment variables or defaults).
    # redo: bool, optional
    #     If True, the function will look for files in the 'sosredo' directory instead of the standard 'sos' directory (default is False).
    # mask_end: bool, optional
    #     If True, the function will mask out the ends of the spectra (outside of the wavelength range of 5000-10000 Angstroms) 
    #       during plotting (default is False).
    # ToOs: bool, optional
    #     If True, only fibers that are marked as ToO (Target of Opportunity) in the fiber mapping will be plotted (default is False).
    # assigned: bool, optional
    #     If True, only fibers that are marked as assigned to targets in the fiber mapping will be plotted (default is False).
    # science: bool, optional
    #     If True, only fibers that are marked as science targets in the fiber mapping will be plotted (default is False).
    # pdf: bool, optional
    #     If True, all selected fibers will be plotted into a single multi-panel PDF instead of individual PNG files (default is True).
    # single_ccd: bool, optional
    #     If True, only the specified CCD will be plotted, and the function will look for science frame files corresponding to that CCD only 
    #       (default is False, in which case it will look for all CCDs in the exposure).
    
    mjd = str(mjd)
    if sos_dir is None:
        sos_dir = os.getenv('BOSS_SOS_S') if obs.lower() == 'lco' else os.getenv('BOSS_SOS_N')
    if sos_dir is None:
        sos_dir = '/data/boss/sos'
        if redo:
            sos_dir = '/data/boss/sosredo'

    outdir = ptt.join(sos_dir,f'{mjd}','plots')

    
    if not single_ccd:
        sci_f = ptt.join(sos_dir,mjd,f'sci-*-??-{str(expid).zfill(8)}.fits')
        if len(glob.glob(sci_f)) == 0:
            splog.warning('No science file found for expid='+str(expid))
            return
        ccds = {x.split('-')[-2][0]: x for x in glob.glob(sci_f) if x.split('-')[-2][0] in ['b','r']}
    else:
        sci_f = ptt.join(sos_dir,mjd,f'sci-*-{ccd_name}?-{str(expid).zfill(8)}.fits')
        if len(glob.glob(sci_f)) == 0:
            splog.warning('No science file found for expid='+str(expid))
            return
        ccds = [ccd_name]

    spectra_dict = {}
    for ccd in ccds:
        spectra_dict[ccd] = SciFrame(mjd, ccd, ccds[ccd], sos_dir)
    
    mask = None
    s = spectra_dict[next(iter(spectra_dict))]
    mask = np.full((s.sci.shape[0],), True)
    if ToOs:
        mask = mask &(s.fmap['TOO'] == 1)
    if assigned:
        mask = mask & (s.fmap['ASSIGNED'] == 1)
    if science:
        mask = mask & (s.fmap['OBJTYPE'] == 'science')

    if mask is None:
        fiberids = np.arange(s.sci.shape[0])
    else:
        fiberids = np.where(mask)[0]

    if len(fiberids) == 0:
        splog.warning('No fibers found for ToO/assigned criteria for expid='+str(expid))
        return
    
    os.makedirs(outdir,exist_ok=True)
    expid = str(expid).zfill(8)

    if pdf:
        num_panels = len(fiberids)
        panels_per_page = 6
        num_pages = (num_panels + panels_per_page - 1) // panels_per_page  # Calculate total number of pages

        # Create a PDF file
        with PdfPages(ptt.join(outdir,f'{expid}-{mjd}_{ccd_name}.pdf')) as pdf:
            fiber = 0
            for page in range(num_pages):
                fig, axes = plt.subplots(6, 1, figsize=(8.5,11)) 
                axes = axes.flatten()  # Flatten to 1D array for easier iteration

                for i in range(panels_per_page):
                    panel_idx = page * panels_per_page + i
                    if panel_idx < num_panels:
                        ax = axes[i]
                        _ = _plotone(ax, spectra_dict, fiberids, fiber, mask_end=mask_end)
                        fiber += 1
                    else:
                        axes[i].axis('off')  # Hide unused subplots
                fig.text(0.01, 0.01, f'ExpID: {int(expid)} MJD: {mjd}', fontsize=8, ha='left', va='bottom')
                fig.text(0.99,0.99, f'Generated on {datetime.now().strftime("%Y-%m-%d %H:%M")} with RUN2D={spectra_dict[ccd].run2d}', fontsize=8, ha='right', va='top')
                plt.tight_layout(pad=1.0, rect=[0.03, 0.03, 0.97, 0.97])
                pdf.savefig(fig)
                plt.close(fig)
        if not single_ccd:
            if lock(ptt.join(outdir,f'{expid}-{mjd}.pdf'), pause=10, niter=4):
                try:
                    os.rename(ptt.join(outdir,f'{expid}-{mjd}_{ccd_name}.pdf'), ptt.join(outdir,f'{expid}-{mjd}.pdf'))
                    splog.info('Saved to '+ptt.join(outdir,f'{expid}-{mjd}.pdf'))
                finally:
                    unlock(ptt.join(outdir,f'{expid}-{mjd}.pdf'))
            else:
                splog.warning('Could not acquire lock for '+ptt.join(outdir,f'{expid}-{mjd}.pdf'))
                splog.info('Saved to '+ptt.join(outdir,f'{expid}-{mjd}_{ccd_name}.pdf'))
        else:
            splog.info('Saved to '+ptt.join(outdir,f'{expid}-{mjd}_{ccd_name}.pdf'))
    else:
        files = Table(names=('name', 'RA', 'DEC', 'title'), dtype=(str, float, float, int))  
        for i, fiber in enumerate(fiberids):
            fig, ax = plt.subplots(figsize=(10,6))
            ra, dec, title = _plotone(ax, spectra_dict, fiberids, i, mask_end=mask_end)
            fig.text(0.01, 0.01, f'ExpID: {int(expid)} MJD: {mjd}', fontsize=8, ha='left', va='bottom')
            fig.text(0.99,0.99, f'Generated on {datetime.now().strftime("%Y-%m-%d %H:%M")} with RUN2D={s.run2d}', fontsize=8, ha='right', va='top')
            plt.tight_layout(pad=1.0, rect=[0.03, 0.03, 0.97, 0.97])
            
            if not single_ccd:
                filename = f'{expid}-{mjd}-{fiber}.png'
    
            else:
                filename = f'{expid}-{mjd}_{ccd_name}-{fiber}.png'

            if lock(ptt.join(outdir,filename), pause=10, niter=4):
                try:
                    plt.savefig(ptt.join(outdir,filename))
                    #splog.info('Saved individual fiber plot to '+ptt.join(outdir,filename))
                finally:
                    unlock(ptt.join(outdir,filename))
            else:
                splog.warning('Could not acquire lock for '+ptt.join(outdir,filename))
                filename = f'{expid}-{mjd}_{ccd_name}-{fiber}.png'
                splog.info('Saved individual fiber plot to '+ptt.join(outdir,filename))
                plt.savefig(ptt.join(outdir,filename))

            plt.close(fig)
            files.add_row((ptt.basename(filename).replace('.png',''), ra, dec, title))
        
        template = ptt.join(idlspec2d_dir, 'templates', 'html','spec_index.html')
        jinja_data = dict(
            pmjd=f'{mjd}-{expid}',
            files = files,
            lsdr10=True, use_thumbs=False)
        

        if lock(ptt.join(outdir, 'tmp-index.html'), pause=10, niter=4):
            try:
                with open(ptt.join(outdir, 'tmp-index.html'), 'w') as f:
                    with open(template) as t:
                        html = Template(t.read()).render(jinja_data)
                        f.write(html)
                os.rename(ptt.join(outdir, 'tmp-index.html'), ptt.join(outdir, f'{expid}-{mjd}.html'))
                splog.info('Saved individual fiber plots to '+ptt.join(outdir, f'{expid}-{mjd}.html'))
            finally:
                unlock(ptt.join(outdir, 'tmp-index.html'))
        else: 
            splog.warning('Could not acquire lock for '+ptt.join(outdir, 'tmp-index.html'))
        

        splog.info('Saved individual fiber plots to '+outdir)
        