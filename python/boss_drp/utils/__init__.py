
from .match import match
from .grep import grep
from .get_lastline import get_lastline
from .load_module import load_module, load_env
from .get_dirs import get_dirs
from .mjd_match import mjd_match
from .splog import Splog
from .merge_dm import merge_dm
from .find_nearest_indx import find_nearest_indx


import astropy.time
import datetime

jdate = str(int(float(astropy.time.Time(datetime.datetime.utcnow()).jd) - 2400000.5))
