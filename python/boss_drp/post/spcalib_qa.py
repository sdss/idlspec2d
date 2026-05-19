from boss_drp import idlspec2d_dir
from boss_drp.utils.splog import splog
from boss_drp.field import Field
from boss_drp.field.field_to_string import field_to_string
from boss_drp.utils.cpbackup import cpbackup
from boss_drp.utils.lock import lock, unlock
from boss_drp.summary import summary_names, Summary_names
from boss_drp.utils import jdate, retry
from boss_drp.utils.parquet.write import write_parquet
from boss_drp.utils.parquet.schema import Schema


from pydl.pydlutils import sdss

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit
import os
import os.path as ptt
from pathlib import Path
from astropy.io import fits
from astropy.table import Table, vstack
from datetime import datetime
import warnings


maskbitfile = ptt.join(os.getenv("IDLUTILS_DIR"),"data","sdss","sdssMaskbits.par")
sdss.set_maskbits(maskbits_file=maskbitfile)

def gauss1(x, amp, mean, sigma):
    return amp * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

def std_hist(fratio, bs=0.015, xmin=-0.3, xmax=-0.2, filt='g', ax = None, run2d=None, lab_loc = (.05,.78), color='k'):
    logf = np.log10(fratio)
    mask = (logf >= xmin) & (logf <= xmax)
    
    if np.count_nonzero(mask) == 0:
        print('test')
        splog.error("catastrophic failure of Flux calibration")
        return np.nan, np.nan

    hist, bin_edges = np.histogram(logf, bins=np.arange(xmin, xmax + bs, bs))
    bin_centers = bin_edges[:-1] + bs / 2

    try:
        p0 = [1, 0.2, 1]
        popt, _ = curve_fit(gauss1, bin_centers, hist, p0=p0,maxfev=10000,
                            bounds=([0, -np.inf, 1e-6], [np.inf, np.inf, np.inf]))
    except RuntimeError:
        try:
            # Use the histogram peak as the initial amplitude and mean
            amp0 = hist.max()
            mean0 = bin_centers[np.argmax(hist)]

            # Rough width estimate from the data near the peak
            sigma0 = np.std(np.repeat(bin_centers, np.maximum(hist.astype(int), 0)))
            if not np.isfinite(sigma0) or sigma0 <= 0:
                sigma0 = 0.05  # fallback

            p0 = [amp0, mean0, sigma0]
            popt, _ = curve_fit(gauss1, bin_centers, hist, p0=p0,maxfev=10000,
                                bounds=([0, -np.inf, 1e-6], [np.inf, np.inf, np.inf]))
        except RuntimeError:
            print('test')
            popt = [np.nan, np.nan, np.nan]

    label = 'Data' if run2d is None else f'Data {run2d}'
    ax.plot(bin_centers, hist, drawstyle='steps-mid', ls='-', lw=1, color=color, label=label)
    ax.set_xlabel(f'log (synflux/calibflux) [{filt}]')
    ax.set_ylabel('Nstd')
    ax.text(lab_loc[0], lab_loc[1], f'mean, sigma\n{popt[1]:.2f}, {popt[2]:.3f}', fontsize=12,
            transform=ax.transAxes, color=color)
    #ax.text(xmin + 0.02, 0.88 * max(hist), f'mean, sigma\n{popt[1]:.2f}, {popt[2]:.3f}', fontsize=12)

    label = 'Fit' if run2d is None else f'Fit {run2d}'
    xfit = np.arange(xmin, xmax, 0.001)
    ax.plot(xfit, gauss1(xfit, *popt), ls='--', color=color, label=label, lw = .5)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(axis='x', which='both', direction='in', top=True)
    ax.tick_params(axis='y', which='both', direction='in', right=True)

    return popt[1], popt[2]  # mean, sigma

def plot_std_hists(f_ratio, ind_good, outname_ps, fieldid=None, mjd=None,
                   f_ratio_alt=None, valid_alt=None, run2d=None, run2d_alt=None):
    if len(ind_good) == 0:
        fit_g = [np.nan, np.nan]
        fit_r = [np.nan, np.nan]
        fit_i = [np.nan, np.nan]
        return fit_g, fit_r, fit_i

    bs = 0.015
    xmin, xmax = -0.3, 0.2

    fig, axs = plt.subplots(3, 1, figsize=(11, 8.5))  # landscape size

    fit_g = std_hist(f_ratio[ind_good,1], bs=bs, xmin=xmin, xmax=xmax, filt='g', ax=axs[0], run2d=run2d)
    fit_r = std_hist(f_ratio[ind_good,2], bs=bs, xmin=xmin, xmax=xmax, filt='r', ax=axs[1], run2d=run2d)
    fit_i = std_hist(f_ratio[ind_good,3], bs=bs, xmin=xmin, xmax=xmax, filt='i', ax=axs[2], run2d=run2d)

    if f_ratio_alt is not None:
        kwrds = dict(run2d=run2d_alt, color='r', lab_loc = (.05,.58))
        _ = std_hist(f_ratio_alt[valid_alt,1], bs=bs, xmin=xmin, xmax=xmax, filt='g', ax=axs[0], **kwrds)
        _ = std_hist(f_ratio_alt[valid_alt,2], bs=bs, xmin=xmin, xmax=xmax, filt='r', ax=axs[1], **kwrds)
        _ = std_hist(f_ratio_alt[valid_alt,3], bs=bs, xmin=xmin, xmax=xmax, filt='i', ax=axs[2], **kwrds)

    if fieldid is not None:
        title = f'Field={fieldid} MJD={str(mjd)}'
        fig.suptitle(title, fontsize=14, y=0.95)

    for ax in axs:
        #ax.grid(True)
        ax.set_xlim(xmin, xmax)
    
    axs[0].legend()

    fig.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(outname_ps, bbox_inches='tight')
    plt.close(fig)

    return fit_g, fit_r, fit_i


def update_spcalib_qa(fieldid, mjd, fit_g, fit_r, fit_i, ct_std, out_fits, spall, run2d, to_fits = True):
    fieldid = int(fieldid)
    obs = spall[0]['OBS'] if 'OBS' in spall.colnames else 'APO'
    mjd = int(mjd)
    if lock(out_fits, pause=10, niter=12):
        try:
            if os.path.exists(out_fits):
                ins = Table.read(out_fits)
                # with fits.open(out_fits, memmap=False) as hdul:
                #     ins = Table(hdul[1].data)
                
                match = np.where((ins['FIELD'] == fieldid) & (ins['MJD'] == mjd))[0]

                new_row = Table(rows=[[fieldid, mjd, obs,
                                       fit_g[0], fit_g[1],
                                       fit_r[0], fit_r[1],
                                       fit_i[0], fit_i[1],
                                       ct_std]],
                                names=('FIELD', 'MJD', 'OBS',
                                       'G_MEAN', 'G_SIG',
                                       'R_MEAN', 'R_SIG',
                                       'I_MEAN', 'I_SIG',
                                       'N_STD'))
                if len(match) == 0:
                    outs = vstack([ins, new_row])
                else:
                    for col in new_row.colnames:
                        ins[col][match[0]] = new_row[col][0]
                    outs = ins
            else:
                outs = Table(rows=[[fieldid, mjd, obs,
                                    fit_g[0], fit_g[1],
                                    fit_r[0], fit_r[1],
                                    fit_i[0], fit_i[1],
                                    ct_std]],
                             names=('FIELD', 'MJD', 'OBS',
                                    'G_MEAN', 'G_SIG',
                                    'R_MEAN', 'R_SIG',
                                    'I_MEAN', 'I_SIG',
                                    'N_STD'))
                
            try:
                date = datetime.now(datetime.timezone.utc).isoformat()
            except:
                date = datetime.utcnow().isoformat()

            schema = Schema()
            schema.yanny_file = ptt.join(idlspec2d_dir, 'datamodel', 'spcalib_qa_dm.par')
            schema.datamodel_arrow_schema('EXT1')
            schema.datamodel_header_metadata('HDR0')
            write_parquet(outs, Path(out_fits.replace('.fits.gz','.parquet')), None, schema=schema, 
                          vo=True, metadata=dict(run2d = run2d, Date = date))

            if to_fits:
                out_fits = out_fits.replace('.parquet','.fits.gz')
                hdu = fits.BinTableHDU(outs, name='spcalib_qa')
                hdu.header.comments['EXTNAME'] = 'SpectroPhotometricQA'
                prim = fits.PrimaryHDU()
                prim.header['RUN2D'] = (run2d, 'IDLSPEC2D RUN2D Version')
                prim.header['DATE'] = (date, 'File creation date (UTC)')

                cols = {'FIELD': 'SDSS FieldID (plateID for plate era data)',
                        'MJD': 'Modified Julian date of combined Spectra',
                        'OBS': 'Observatory of Observation',
                        'G_MEAN': 'mean(synflux/calibflux) [g]',
                        'G_SIG':  'sigma(synflux/calibflux) [g]',
                        'R_MEAN': 'mean(synflux/calibflux) [r]',
                        'R_SIG':  'sigma(synflux/calibflux) [r]',
                        'I_MEAN': 'mean(synflux/calibflux) [i]',
                        'I_SIG':  'sigma(synflux/calibflux) [i]',
                        'N_STD':  'Number of Spectro-Photometric Standards'}
                
                for card in hdu.header.cards:
                    if str(card[1]).strip().upper() in cols.keys():
                        hdu.header[card[0]] = (card[1], cols[card[1]])
                
                hdulist = fits.HDUList([prim, hdu])

                
                hdulist.writeto(out_fits, overwrite=True)

        finally:
            unlock(out_fits)
    return



def spcalib_qa(run2d=None, field=None, mjd=None, rerun=False,
                catchup=False, nobkup=False, epoch=False, 
                boss_spectro_redux = None, outdir = None, to_fits=True, 
                run2d_alt=None, boss_spectro_redux_alt=None, plot_only=False):
    log_folder = '.'
    splog.open(ptt.join(log_folder, jdate.astype(str)+'.log'))

    run2d = run2d or os.getenv('RUN2D')

    flag = '-epoch' if epoch else ''

    outname = f'spCalib_QA-{run2d}{flag}'
    if run2d_alt is not None:
        outname = f'spCalib_QA-{run2d}-{run2d_alt}{flag}'

    if boss_spectro_redux is not None:
        os.environ['BOSS_SPECTRO_REDUX'] = boss_spectro_redux
    summary_names.build(os.getenv('BOSS_SPECTRO_REDUX'), run2d=run2d, epoch=epoch)
    boss_root = summary_names.outdir

    out_fits = ptt.join(boss_root, f'{outname}.fits')
    spall_full_file = ptt.join(boss_root, f'spAll-{run2d}{flag}.parquet')
#    spall_full_file = ptt.join(boss_root, f'spAll-{run2d}{flag}.fits.gz')

    if run2d_alt is not None:
        if boss_spectro_redux_alt is None:
            boss_spectro_redux_alt = os.getenv('BOSS_SPECTRO_REDUX')
        sna = Summary_names()
        sna.build(boss_spectro_redux_alt, run2d=run2d_alt, epoch=epoch)
        boss_root_a = sna.outdir
        spall_full_file_a = ptt.join(boss_root_a, f'spAll-{run2d_alt}{flag}.parquet')
        #spall_full_file_a = ptt.join(boss_root_a, f'spAll-{run2d_alt}{flag}.fits.gz')

    if field is not None:
        spallfile=f'spAll-{field_to_string(field)}-{mjd}.fits'
        fs = Field(os.getenv('BOSS_SPECTRO_REDUX'),run2d,field,mjd=mjd,epoch=epoch)
        spallfile = ptt.join(fs.spec_dir(), spallfile)

        if run2d_alt:
            spallfile_a=f'spAll-{field_to_string(field)}-{mjd}.fits'
            fs_a = Field(boss_spectro_redux_alt,run2d_alt,field,mjd=mjd,epoch=epoch)
            spallfile_alt = ptt.join(fs_a.spec_dir(), spallfile_a)
        
        outname = f'{outname}-{field_to_string(field)}-{mjd}'
        outname_ps = ptt.join(fs.dir(),outname+'.png')
        logfile = ptt.join(fs.dir(),outname+'.log')
        if outdir:
            outname_ps = ptt.join(outdir,outname+'.png')
            logfile = ptt.join(outdir,outname+'.log')

        
        
        splog.open(logfile=logfile, backup=(not nobkup))
        
        splog.info(f'Log file {logfile} opened {datetime.now().strftime("%a %b %d %H:%M:%S %Y")}')

        if not nobkup: cpbackup(outname_ps)
    else:
        if not rerun:
            if catchup:
                # Load full spAll FITS table
                if not ptt.exists(spall_full_file):
                    if not ptt.exists(spall_full_file.replace('.parquet','.fits.gz')):
                        splog.warning(f"Missing spAll file: {spall_full_file}")
                        return
                    else:
                        spall_full_file = spall_full_file.replace('.parquet','.fits.gz')
                spall_tag = Table.read(spall_full_file)

                fieldids = spall_tag['FIELD'].data
                mjds = spall_tag['MJD'].data
                
                _, idx = np.unique(fieldids, return_index=True)
                ufieldids = np.sort(idx)  # optional, if you want them in order of appearance
                unique_fieldids = fieldids[ufieldids]
                
                if ptt.exists(out_fits): 
                    ins_e = Table.read(out_fits)
                    #ins_e = fits.getdata(out_fits, 1)
                for i in range(len(unique_fieldids)):
                    current_field = unique_fieldids[i]
                    ids = np.where(fieldids == current_field)[0]
                    
                    if len(ids) > 0:
                        matching_mjds = mjds[ids]
                        _, idx_mjd = np.unique(matching_mjds, return_index=True)
                        unique_mjds = matching_mjds[np.sort(idx_mjd)]

                        for j in range(len(unique_mjds)):
                            # Optional: check against existing entries if `ins_e` is set
                            if ins_e is not None:
                                match = np.where((ins_e['FIELD'] == current_field) & (ins_e['MJD'] == unique_mjds[j]))[0]
                                if len(match) > 0:
                                    splog.info(f"Skipping Existing: {str(current_field).strip()}-{str(unique_mjds[j]).strip()}")
                                    continue

                            # Call processing routine
                            spcalib_qa(run2d=run2d, field=current_field,
                                       mjd=unique_mjds[j], nobkup=nobkup,
                                       epoch=epoch, run2d_alt=run2d_alt,plot_only=plot_only,
                                       boss_spectro_redux_alt=boss_spectro_redux_alt)
            outname_ps = ptt.join(boss_root, outname_ps)

        else:
            # Load full spAll FITS table
            if not ptt.exists(spall_full_file):
                if not ptt.exists(spall_full_file.replace('.parquet','.fits.gz')):
                    splog.warning(f"Missing spAll file: {spall_full_file}")
                    return
                else:
                    spall_full_file = spall_full_file.replace('.parquet','.fits.gz')
            spall_tag = Table.read(spall_full_file) #fits.getdata(spall_full_file)
            fieldids = spall_tag['FIELD']
            mjds = spall_tag['MJD']
            
            # Get unique field IDs in order of first appearance
            _, idx_field = np.unique(fieldids, return_index=True)
            ufieldids = np.sort(idx_field)
            unique_fieldids = fieldids[ufieldids]

            # Loop through unique fields and their unique MJDs
            for i in range(len(unique_fieldids)):
                current_field = unique_fieldids[i]
                ids = np.where(fieldids == current_field)[0]

                if len(ids) > 0:
                    matching_mjds = mjds[ids]
                    _, idx_mjd = np.unique(matching_mjds, return_index=True)
                    unique_mjds = matching_mjds[np.sort(idx_mjd)]

                    for j in range(len(unique_mjds)):
                        spcalib_qa(run2d=run2d, field=current_field,
                                   mjd=unique_mjds[j], nobkup=nobkup, 
                                   epoch=epoch,run2d_alt=run2d_alt, plot_only=plot_only,
                                   boss_spectro_redux_alt=boss_spectro_redux_alt)

        
        spallfile = spall_full_file
        spallfile_alt = spall_full_file_a
        outname_ps = ptt.join(boss_root, outname+'.png')
        
    # Load full spAll FITS table
    if not ptt.exists(spallfile):
        if ptt.exists(spallfile+'.gz'):
            spallfile = spallfile+'.gz'
        else:
            splog.warning(f"Missing spAll file: {spallfile}")
            return

    spall_data = retry(Table.read, spallfile, retries = 3, delay = 5, noerr=True,
                        logger=splog.info)
    # spall_data = retry(fits.getdata, retries = 3, delay = 5, noerr=True,
    #                     logger=splog.info, filename = spallfile)
    if spall_data is None:
        splog.warning('Skipping {spallfile}')
        return

    # Filter by standards
    is_std = (spall_data['OBJTYPE'] == 'SPECTROPHOTO_STD') & ((spall_data['ZWARNING'] &
              sdss.sdss_flagval('ZWARNING',['UNPLUGGED']))==0)
    
    ct_std = np.count_nonzero(is_std)
    if ct_std == 0:
        splog.error(f"No standards found in {spallfile}")
        fit_g = [np.nan, np.nan]
        fit_r = [np.nan, np.nan]
        fit_i = [np.nan, np.nan]
    else:
        splog.info(f"Loaded {np.count_nonzero(is_std)} standard stars from {spallfile}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            f_ratio = spall_data['SPECTROSYNFLUX'][is_std] / spall_data['CALIBFLUX'][is_std]

        # Only use values with valid flux ratio
        valid = f_ratio[:,1] > 0

        if run2d_alt is None:
            fit_g, fit_r, fit_i = plot_std_hists(f_ratio, valid, outname_ps,
                                                 fieldid=field, mjd=mjd)
    

    if run2d_alt is not None:
        if not ptt.exists(spallfile_alt):
            if ptt.exists(spallfile_alt+'.gz'):
                spallfile_alt = spallfile_alt+'.gz'
            else:
                splog.warning(f"Missing spAll file: {spallfile_alt}")
                return

        spall_data_alt = retry(Table.read, spallfile_alt, retries = 3, delay = 5, noerr=True,
                        logger=splog.info)
        # spall_data_alt = retry(fits.getdata, retries = 3, delay = 5, noerr=True,
        #                 logger=splog.info, filename = spallfile_alt)
        if spall_data_alt is None:
            splog.warning('Skipping {spallfile_alt}')
            return

        # Filter by standards
        is_std_a = (spall_data_alt['OBJTYPE'] == 'SPECTROPHOTO_STD') & ((spall_data_alt['ZWARNING'] &
                sdss.sdss_flagval('ZWARNING',['UNPLUGGED']))==0)
        
        ct_std_a = np.count_nonzero(is_std_a)
        if ct_std_a == 0:
            splog.error(f"No standards found in {spallfile_alt}")
            fit_g_a = [np.nan, np.nan]
            fit_r_a = [np.nan, np.nan]
            fit_i_a = [np.nan, np.nan]
        else:
            splog.info(f"Loaded {np.count_nonzero(is_std_a)} standard stars from {spallfile_alt}")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                f_ratio_a = spall_data_alt['SPECTROSYNFLUX'][is_std_a] / spall_data_alt['CALIBFLUX'][is_std_a]

            # Only use values with valid flux ratio
            valid_a = f_ratio_a[:,1] > 0

        fit_g, fit_r, fit_i = plot_std_hists(f_ratio, valid, outname_ps,
                                             fieldid=field, mjd=mjd,
                                             f_ratio_alt = f_ratio_a, 
                                             valid_alt = valid_a,
                                             run2d=run2d, run2d_alt = run2d_alt)

    if not plot_only:
        if field is not None:
            update_spcalib_qa(field, mjd, fit_g, fit_r, fit_i, len(valid), out_fits, 
                              spall_data, run2d, to_fits=to_fits)

    splog.info('SpectroPhoto QA Complete')
    splog.close()


if __name__ == '__main__':
    spcalib_qa(run2d='v6_2_1-fix', field=104648, mjd=59984, rerun=False,
                catchup=False, nobkup=False, epoch=False)
