from .transform import transform
from .html import html
from .plan import Plan
from .report_err import Report

import matplotlib
import os
import os.path as ptt
import time

def boss_arcs_to_traces(mjd = None, outdir = None, obs = 'lco', vers = 'master',
                        threads = 8, nskip = 40, cams = None, fitsname = None,
                        designMode = 'uknown', sosdir = None, clobber = False):
    if vers.lower() == 'sos':
        vers = ''
        sos = True
        if sosdir is not None:
            os.environ["BOSS_SPECTRO_REDUX"]  = ptt.join(f'{sosdir}',f'{mjd}')
    else:
        sos = False
        
    try:
        import logging
        _log = logging.getLogger("astropy")
        _log.setLevel(logging.CRITICAL)
        mjd = int(mjd)
        transform(mjd, obs=obs, clobber=clobber, threads=threads,
                 outdir=outdir, vers=vers, cams=cams, sos=sos,
                 designMode=designMode)
    except Exception as e:
        import traceback
        print(traceback.format_exc())
        print(type(e).__name__, ":", e)
        if not sos:
            exit()

    if sos:
        soshtml(mjd, obs, sosdir)

    else:
        print(f'Successful completion of boss_arcs_to_trace for mjd {mjd} obs {obs} at {time.ctime()}')
