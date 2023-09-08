
import h5py
import healpy as hp
import numpy as np
import json
from os import getenv
import os.path as ptt

def find_nearest_idx(value,array):
    value = np.atleast_1d(value)
    array = np.asarray(array)
    idx = np.full(len(value), -1)
    for i, val in enumerate(value):
        idx[i] = (np.abs(array - val)).argmin()
    return idx

class simple_dust_2023:
    def __init__(self, mapfile = None):
        self.load_map(mapfile=mapfile)
    
    def load_map(self,mapfile=None):
        if mapfile is not None:
            if not ptt.exists(mapfile):
                print(f'{mapfile} not found using default')
        if mapfile is None:
            with open(getenv('DUSTMAPS_CONFIG_FNAME')) as config:
                conf = json.load(config)
            mapfile = ptt.join(conf['data_dir'],'manual','map3d_multires.h5')
        f1 = h5py.File(mapfile,'r')
        self._distance_centers = f1['distance_centers'][:]
        self._extinction = f1['extinction'][:]
        self._extinction_variance = f1['extinction_variance'][:]
        self._nside = f1['nside'][:]

    def _find_data_idx(self, b,l,nside=512,nest=True):
        if type(l) is list:
            l = np.asarray(l)
        if type(b) is list:
            b = np.asarray(b)
        theta = np.radians(90. - b)
        phi = np.radians(l)
        
        if not hasattr(l, '__len__'):
            if (b < -90.) or (b > 90.):
                return -1

            pix_idx = hp.pixelfunc.ang2pix(nside, theta, phi, nest=nest)

            return pix_idx

        idx = (b >= -90.) & (b <= 90.)

        pix_idx = np.empty(l.shape, dtype='i8')
        pix_idx[idx] = hp.pixelfunc.ang2pix(nside, theta[idx], phi[idx], nest=nest)
        pix_idx[~idx] = -1
        return pix_idx



    def query(self, coords):
        # Get number of coordinates requested
        n_coords_ret = coords.shape[0]

        # Determine if distance has been requested
        has_dist = hasattr(coords.distance, 'kpc')
        d = coords.distance.kpc

        # Extract the correct angular pixel(s)
        # t0 = time.time()
        pix_idx = self._find_data_idx(coords.l.deg, coords.b.deg)
        in_bounds_idx = np.where((pix_idx != -1))[0]
        in_dist_idx = np.where((d > 0.0) & (d < self._distance_centers[-1]+(self._distance_centers[-1]-self._distance_centers[-2])/2))[0]

        ret = np.full((n_coords_ret,), np.nan, dtype='f4')
        ret = self._extinction[find_nearest_idx(d, self._distance_centers),pix_idx]
        ret[~in_bounds_idx] = np.nan
        ret[~in_dist_idx] = np.nan

        return ret

