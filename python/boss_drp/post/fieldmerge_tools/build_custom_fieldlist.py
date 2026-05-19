from boss_drp.field import field_to_string, Field, fieldgroup
from boss_drp.utils import  get_lastline

from glob import glob
import os.path as ptt
from astropy.table import Table, vstack
from astropy.io import fits
import numpy as np
import warnings

def build_custom_fieldlist(indir, custom, run2d, run1d):
    flist = Table()
    spfields = []
    for cc in [custom,custom+'_lco',custom+'_apo']:
        fc = Field(indir, run2d, cc, custom_name = cc, custom = True)
        spfields.extend(glob(ptt.join(fc.dir(), f'spFullsky-{cc}-*.fits')))
    for spfield in spfields:
        hdr = fits.getheader(spfield,0)
        
        cc = ptt.basename(spfield).split('-')[1]
        if cc.split('_')[-1] in ['lco','apo']:
            obs = cc.split('_')[-1].upper()
        else:
            obs = ''
        fc = Field(indir, run2d, cc, custom_name = cc, custom = True)
        spdiagcomblog1 = ptt.join(fc.dir(),f"spDiagcomb-{cc}-{str(hdr['RUNMJD'])}.log")
        spdiagcomblog = ptt.join(fc.dir(),f"spDiagcomb-{cc}-{str(hdr['RUNMJD'])}_{str(hdr['MJD'])}.log")

        if ptt.exists(spdiagcomblog):
            lastline = get_lastline(spdiagcomblog)
            if 'Successful completion' in lastline:
                #Case where this 1D log file completed, which is not a case that should ever occur
                STATUSCOMBINE = 'Done'
            else:
                #Case where this 1D log file isn't completed
                STATUSCOMBINE = 'RUNNING'
        elif ptt.exists(spdiagcomblog1):
            lastline = get_lastline(spdiagcomblog1)
            if 'Successful completion' in lastline:
                #Case where this 1D log file completed, which is not a case that should ever occur
                STATUSCOMBINE = 'Done'
            else:
                #Case where this 1D log file isn't completed
                STATUSCOMBINE = 'RUNNING'
        else:
            STATUSCOMBINE = 'Pending'
        spDiag1dlog = ptt.join(fc.spec1d_dir(run1d), f"spDiag1d-{cc}-{str(hdr['MJD'])}.log")
        if ptt.exists(spDiag1dlog):
            lastline = get_lastline(spDiag1dlog)
            if 'Successful completion' in lastline:
                #Case where this 1D log file completed, which is not a case that should ever occur
                STATUS1D = 'Done'
            else:
                #Case where this 1D log file isn't completed
                STATUS1D = 'RUNNING'
        else:
            STATUS1D = 'Pending'

        try:
            sn2_g1 = hdr['SPEC1_G']
        except:
            sn2_g1 = np.NaN
        try:
            sn2_i1 = hdr['SPEC1_I']
        except:
            sn2_i1 = np.NaN

        try:
            sn2_g2 = hdr['SPEC2_G']
        except:
            sn2_g2 = np.NaN
        try:
            sn2_i2 = hdr['SPEC2_I']
        except:
            sn2_i2 = np.NaN
            
        if obs.lower() == 'lco':
            sn2 = [sn2_g2,sn2_i2]
        elif obs.lower() == 'apo':
            sn2 = [sn2_g1,sn2_i1]
        else:
            sn2 = [sn2_g1,sn2_i1,sn2_g2,sn2_i2]
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')

            flist_row = Table({'RUN2D':[run2d], 'RUN1D':[run1d],
                               'PROGRAMNAME':[custom], 'FIELD':[0], 'MJD': [hdr['MJD']],
                               'STATUS2D':['Done'], 'STATUSCOMBINE':[STATUSCOMBINE],
                               'STATUS1D':[STATUS1D],'FIELDQUALITY':['good'],
                               'FIELDSN2':[np.nanmin( np.array(sn2))],
                               'OBSERVATORY':[obs]})

        flist =  vstack([flist, flist_row])
    flist.pprint_all()
    return(flist)

