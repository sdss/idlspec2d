from astropy.coordinates import SkyCoord
import astropy.units as units
from dustmaps.bayestar import BayestarWebQuery
import numpy as np
import sys

sys.argv=filter(None,sys.argv) 
ll=np.float(sys.argv[1])
bb=np.float(sys.argv[2])
rr=np.float(sys.argv[3])
bayestar = BayestarWebQuery(version='bayestar2017')
coords = SkyCoord(90.*units.deg, 30.*units.deg,distance=100.*units.pc, frame='galactic')
reddening = bayestar(coords, mode='median')
print reddening