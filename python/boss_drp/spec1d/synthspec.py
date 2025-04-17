from boss_drp.utils.splog import splog

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

def synthspec(zans, loglam, eigendir=None):
    """
    Generate synthetic spectra from SDSS eigentemplates for one or more objects.
    
    This function can handle single row inputs, lists/arrays of rows, or full
    Astropy Tables. It will generate synthetic spectra based on the redshift fitting
    results contained in the 'zans' argument.

    This is a python port of the idl function of the same name

    Parameters
    ----------
    zans : astropy.table.Row, list of astropy.table.Row, or astropy.table.Table
        A single row or multiple rows (from an Astropy Table) representing redshift
        fitting results (must contain fields like 'TFILE', 'Z', 'TCOLUMN', 'THETA', etc.).
        This function can also accept an entire Astropy Table, in which case it will process
        each row as an individual object.

    loglam : np.ndarray
        Array of log10(wavelength) values to which the spectrum should be interpolated.

    eigendir : str, optional
        Local path to the eigentemplates directory. If not provided, the templates will be downloaded
        from the SDSS GitHub repository.

    Returns
    -------
    newflux : np.ndarray
        The synthetic flux interpolated to the input `loglam` grid, accounting for redshift.
        If no template is available, returns an array of zeros with the same length as `loglam`.

    Notes
    -----
    - If the template file (`TFILE`) is not specified or empty, the function returns a zero flux array.
    - Polynomial terms up to degree `NPOLY` are optionally included in the spectrum model.
    - The final model is redshifted using the provided `Z` value and interpolated to the desired loglam grid.

    """
    
    # Handle multiple zans input
    if isinstance(zans, (list, tuple, np.ndarray, Table)) and not isinstance(zans, Row):
        if isinstance(zans, Table):
            zans_list = zans  # Tables support indexing like lists
        else:
            zans_list = list(zans)
        nobj = len(zans_list)
        out = []
        for iobj in range(nobj):
            print(f'Generating synthetic spectrum for object {iobj + 1}/{nobj}')
            thisloglam = loglam if loglam is None or loglam.ndim == 1 else loglam[:, iobj]
            flux = synthspec(zans_list[iobj], thisloglam, run2d, eigendir=eigendir)
            out.append(flux)
        return np.array(out).T  # shape: (n_pix, n_obj)
    
    if eigendir is None:
        # Try using the local IDLSPEC2D to get templates
        eigendir = os.getenv('IDLSPEC2D_DIR', default=None)
        if eigendir is not None:
            eigendir = os.path.join(eigendir,'templates')

    tfile = str(zans['TFILE']).strip()
    if tfile == '':
        print('No Template')
        return np.zeros(len(loglam), dtype=np.float32)
        
    if eigendir is not None:
        # Use local templates
        model = os.path.join(eigendir, tfile)
    else:
        # Download requested template file
        try:
            baseurl = f'https://raw.githubusercontent.com/sdss/idlspec2d/refs/tags/{run1d}/templates/{tfile}'
            model = download_file(baseurl, cache=True, pkgname='boss_drp')
        except:
            try:
                baseurl = f'https://raw.githubusercontent.com/sdss/idlspec2d/refs/heads/{run1d}/templates/{tfile}'
                model = download_file(baseurl, cache=True, pkgname='boss_drp')
            except:
                baseurl = f'https://raw.githubusercontent.com/sdss/idlspec2d/refs/heads/master/templates/{tfile}'
                model = download_file(baseurl, cache=True, pkgname='boss_drp')  

    try:
        with fits.open(model, ignore_missing_end=True) as hdul:
            # Read the Template file
            starflux = hdul[0].data
            shdr = hdul[0].header
            has_data = starflux is not None and starflux.size > 1
            if not has_data:
                print('Invalid Template File')
                return np.zeros(len(loglam), dtype=np.float32)

            try:
                starloglam0 = shdr['COEFF0']
                stardloglam = shdr['COEFF1']
            except:
                print('Missing Wavelenth Coeffients in Template file')
                return np.zeros(len(loglam), dtype=np.float32)
    except Exception as e:
        print('Error reading Template')
        print(f'{type(e).__name__}: {e}')
        return np.zeros(len(loglam), dtype=np.float32)

    # Reshape and select relevant columns
    starflux = starflux.T if starflux.ndim == 2 else starflux[:, np.newaxis]
    starflux = starflux[:, np.array(zans['TCOLUMN'])[np.where(zans['TCOLUMN'] != -1)[0]]]

    # Add polynomial terms if required
    if zans['NPOLY'] != 0:
        x = np.arange(starflux.shape[0], dtype=np.float32) / starflux.shape[0]
        pa =  np.vstack([x**i for i in range(zans['NPOLY'])]).T  # shape (npix, npoly)
        starflux = np.hstack([starflux, pa])

    # Apply coefficients to get synthetic flux
    synflux = starflux @ zans['THETA'][:starflux.shape[1]]
    
    # Build loglam grid for template
    starloglam = starloglam0 + np.arange(starflux.shape[0]) * stardloglam

    # Redshift and interpolate
    newflux = np.interp(loglam, starloglam + np.log10(1 + zans['Z']), synflux)

    return newflux
