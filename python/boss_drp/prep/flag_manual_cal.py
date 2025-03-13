from boss_drp.field import field_to_string
from boss_drp.utils import jdate

from pydl.pydlutils import yanny
from astropy.table import Table, unique, MaskedColumn
import os.path as ptt
from os import remove, getenv,environ
from collections import OrderedDict

try:
    import git
    if getenv('GIT_PYTHON_TRACE') is None:
        environ['GIT_PYTHON_TRACE'] = '2'  # 'full' gives even more details

except:
    git = None


def flag_manual_cal(type='arc', field=None, mjd=None, obs = 'lco', expid = None, nogit=False):
    spManCal_file = ptt.join(getenv('SDHDRFIX_DIR'), obs.lower(), 'sdHdrfix','spManCal.par')
    
    try:
        spManCal = yanny.read_table_yanny(spManCal_file,'OPMANCAL')
    except :
        spManCal = Table(names=('type','field','mjd','expid'),
                         dtype=('S4','S7',int,int),
                         masked=True,  # Enable masking
                         descriptions=("arc or flat","FieldID", "MJD", "Manual Expid"))
    spManCal.meta =  OrderedDict({'MJD':  jdate.astype(str)+"   # Modified Julian Date of latest edit",
                                  'OBS':  str(obs)       +"     # Observatory" })
    
    if expid is None:
        expid = -999
    else:
        expid = int(expid)
    spManCal.add_row((type.lower(), field_to_string(field), int(mjd), expid))
    spManCal = unique(spManCal, keys=['type','field','mjd'], keep='last')
    spManCal.sort('mjd')
    print('Writing to: '+spManCal_file)
    print(spManCal)
    
    if ptt.exists(spManCal_file):
        remove(spManCal_file)
    yanny.write_ndarray_to_yanny(spManCal_file, spManCal, structnames='OPMANCAL',hdr=spManCal.meta, comments=None)

    if (git is not None) and (not nogit):
        repo = git.Repo(getenv('SDHDRFIX_DIR'))  # Absolute path to repo
        relative_path = ptt.relpath(spManCal_file, getenv('SDHDRFIX_DIR'))  # Convert to relative path
        repo.index.add([relative_path])
        print(f'Adding file to git repo')
    else:
        print('File not yet added to git repo... run the following commands to add it')
        print(f'cd {ptt.dirname(spManCal_file)}')
        print(f'git add {ptt.basename(spManCal_file)}')

def check_manual_cal(type='arc', field=None, mjd=None, obs = 'lco'):
    spManCal_file = ptt.join(getenv('SDHDRFIX_DIR'), obs.lower(), 'sdHdrfix','spManCal.par')
    if not ptt.exists(spManCal_file):
        return False, None
    
    spManCal = yanny.read_table_yanny(spManCal_file,'OPMANCAL')
    spManCal = spManCal[(spManCal['mjd'] == int(mjd)) &
                        (spManCal['field'] == field_to_string(field)) &
                        (spManCal['type'] == type.lower())]
    if len(spManCal) > 0:
        return (True, spManCal['expid'][0] if spManCal['expid'][0] != -999 else None)
    return (False, None)
