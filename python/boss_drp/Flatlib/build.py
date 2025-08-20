
import glob
import os.path as ptt
from os import getenv
from astropy.io import fits
from astropy.coordinates import AltAz, EarthLocation
from astropy.time import Time
import astropy.units as u
from astropy.table import Table
from datetime import datetime, timedelta
from tqdm import tqdm
import numpy as np
from collections import OrderedDict
def build(dir_, obs):
    cols = {'MJD':'MJD','ROT':'ROTPOS','ALT':'ALT','AZ':'AZ','EXP':'exposure','TAI':'TAI','CARTID':'CARTID'}
    catfile = ptt.join(dir_,'calibs','allflats.fits')
    if ptt.exists(catfile):
        catalog = Table.read(catfile, format='fits')
    else:
        catalog = None
    for obs in obs:
        sp= 1 if obs =='apo' else 2
        if obs == 'apo':
            sp = 1
            location = EarthLocation(lat = 32.780361*u.deg, lon=254.179532*u.deg, height=2788*u.m)
            tf = '%Y-%m-%dT%H:%M:%S'
        else:
            location = EarthLocation(lat = -29.01597*u.deg, lon=289.307920*u.deg, height=2380*u.m)
            sp = 2
            tf = '%Y-%m-%dT%H:%M:%S.%f'
        flatlib = ptt.join(dir_,'calibs',obs,'?????',f'spFlat-b{sp}*.fits*')
        for ff in tqdm(glob.glob(flatlib)):
            if catalog is not None:
                if ptt.basename(ff).split('.')[0].split('-')[-1] in catalog['EXP']:
                    continue
            hdr = fits.getheader(ff)
            meta = OrderedDict({})
            for col in cols.keys():
                if col== 'TAI':
                    meta['TAI'] = hdr['TAI-BEG'] + (hdr['EXPTIME']/2.0)
                else:
                    meta[col] = hdr[cols[col]]
                if meta[col] == '':
                    meta[col] = np.NaN
            meta['OBS'] = obs
            meta['qbad_b'] = 1
            meta['qbad_r'] = 1
            obstime = Time(datetime.strptime(hdr['DATE-OBS'], tf) + timedelta(seconds=hdr['EXPTIME']/2.0))
            try:
                altaz = AltAz(alt= meta['ALT']*u.deg,az=meta['AZ']*u.deg,obstime=obstime,location=location)
            except:
                continue
            meta['AIRMASS'] = altaz.secz.value
            if catalog is None:
                for key in meta.keys():
                    meta[key] = [meta[key]]
                catalog = Table(meta)
            else:
                catalog.add_row(meta)
    catalog.write(catfile, overwrite=True)
