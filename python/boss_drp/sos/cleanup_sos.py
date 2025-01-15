from boss_drp.sos.sos_classes import SOS_config
from boss_drp.utils.lock import unlock


import atexit
import signal
import sys
import os.path as ptt
from os import getenv, lstat
from glob import glob
import warnings
import time
import re

def exp_file(filename):
    # Regular expression to match '*-????????.*'
    pattern = r'^.*-\d{8}\..*$'
    return bool(re.match(pattern, filename))

def cleanup():
    files = glob(ptt.join(SOS_config.sosdir,str(SOS_config.MJD),f'*{SOS_config.CCD}*.lock'))
    files.extend(glob(ptt.join(SOS_config.sosdir,str(SOS_config.MJD),'trace',str(SOS_config.MJD),f'*{SOS_config.CCD}*.lock')))
    for f in files:
        if SOS_config.exp is not None:
            if exp_file(ptt.basename(f)):
                if SOS_config.exp not in f:
                    continue
        unlock(f)

    # not checking for /home/sdss5/software/sdsscore/main/*/sdHdrfix/*.lock
    # not checking for any combined (red+blue files) - dont want to unlock if only 1 process is killed for systemd
    
    check()

class FileLockWarning(Warning):
    pass


def check(force_unlock = False):
    print(ptt.join(SOS_config.sosdir,str(SOS_config.MJD),f'*.lock'))
    print(ptt.join(SOS_config.sosdir,str(SOS_config.MJD),'trace',str(SOS_config.MJD),f'*.lock'))
    print(ptt.join(SOS_config.sosdir,'combined',f'*.lock'))
    print(ptt.join(getenv('SDHDRFIX_DIR'),getenv('OBSERVATORY','APO').lower(),'sdHdrfix/*.lock'))
    files = glob(ptt.join(SOS_config.sosdir,str(SOS_config.MJD),f'*.lock'))
    files.extend(glob(ptt.join(SOS_config.sosdir,str(SOS_config.MJD),'trace',str(SOS_config.MJD),f'*.lock')))
    files.extend(glob(ptt.join(SOS_config.sosdir,'combined',f'*.lock')))
    files.extend(glob(ptt.join(getenv('SDHDRFIX_DIR'),getenv('OBSERVATORY','APO').lower(),'sdHdrfix/*.lock')))
    current_time = time.time()
    for file in files:
        if ((current_time - lstat(file).st_ctime) > 300) or (force_unlock):
            warnings.warn(f'Unlocking Locked File: {file}', FileLockWarning)
            unlock(file)
            continue            
        warnings.warn(f'Locked File Exists: {file}',FileLockWarning)
    

# Function to set up the cleanup mechanism
def setup_cleanup():
    # Register cleanup with atexit (runs when the program exits normally)
    atexit.register(cleanup)

    # Handle signals like SIGINT (Ctrl+C) or SIGTERM (kill)
    def signal_handler(sig, frame):
        print(f"Signal {sig} received. Running cleanup.")
        cleanup()
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)  # Handle Ctrl+C
    signal.signal(signal.SIGTERM, signal_handler)  # Handle kill

# Automatically set up cleanup when the module is imported
setup_cleanup()

