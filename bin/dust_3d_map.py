#!/usr/bin/env python
import warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import SkyCoord
import astropy.units as units
from dustmaps.bayestar import BayestarWebQuery
import numpy as np
import sys


sys.argv=filter(None,sys.argv)
ll=np.float(sys.argv[1])
bb=np.float(sys.argv[2])
rr=np.float(sys.argv[3])
bayestar = BayestarWebQuery(version='bayestar2015')
coords = SkyCoord(ll*units.deg, bb*units.deg,distance=rr*units.pc, frame='galactic')
reddening = bayestar(coords, mode='median')#'percentile',pct=70.)#mode='median')
#reddening1 = bayestar(coords, mode='percentile',pct=70.)
#reddening2 = bayestar(coords, mode='percentile',pct=30.)
#reddening=(reddening1+reddening2)/2.0
#reddening=reddening/3.1
if rr <= 00.0:
    reddening=np.nan
print reddening
