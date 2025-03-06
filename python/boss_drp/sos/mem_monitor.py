import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..","..")))
#from utils.splog import splog

from boss_drp.utils.splog import splog
splog._log.setLevel('CRITICAL')

import logging
_log = logging.getLogger("matplotlib.pyplot")
_log.setLevel(logging.CRITICAL)

from matplotlib.pyplot import close as pltclose
from astropy.io.fits import HDUList as fitsHDUList
from astropy.io.fits import FITS_record, FITS_rec

import inspect
import tracemalloc
try:
    import psutil
except:
    psutil = None
import gc
import io
splog._log.setLevel('DEBUG')

def print_usage(deep = False):
    if not tracemalloc.is_tracing():
        return
    current, peak = tracemalloc.get_traced_memory()

    # Get the current call stack
    frame = inspect.currentframe()

    # Find the caller's frame (one level up)
    caller_frame = frame.f_back
    if deep:
        caller_frame = caller_frame.f_back
    try:
        filename = caller_frame.f_code.co_filename
        lineno = caller_frame.f_lineno
    except:
        try:
            caller_frame = frame.f_back
            filename = caller_frame.f_code.co_filename
            lineno = caller_frame.f_lineno
        except:
            return
    splog.debug(f'Memory at {filename}:{lineno} - Current:{current / 10**6} MB - Peak:{peak / 10**6} MB')
    return

def check(prv=None, usage=False):
    if not tracemalloc.is_tracing():
        return [None, None]
    if prv is None: prv = [None,None]
    obj = get_large_objects_gc(threshold=1e6, prv = prv[0])
    files = list_open_files(prv=prv[1])
    if usage:
        print_usage(deep=True)
    return [obj, files]

def get_large_objects_gc(threshold=1e6, prv = None):
    temp = {id(obj): sys.getsizeof(obj) for obj in gc.get_objects() if sys.getsizeof(obj) > threshold}
    if prv is not None:
        leaked_objects = {k: temp[k] for k in temp if k not in prv}
        if leaked_objects:
            for obj_id, obj in leaked_objects.items():
                splog.debug(f"Leaked objects (size in bytes): {obj_id}: {obj} ({type(obj)})")
                for name, ref in globals().items():
                    if id(ref) == obj_id:
                        splog.debug(f"name: {name}")
                        break

        else:
            splog.debug("No large objects leaked.")
        for i,stat in enumerate(tracemalloc.take_snapshot().statistics("lineno")[:10]):
            splog.debug(f'Mem {i}: {stat}')
    return {id(obj): sys.getsizeof(obj) for obj in gc.get_objects() if sys.getsizeof(obj) > threshold}


def list_open_files(prv =None):
    """Returns a list of open file paths for the current process."""
    if psutil is None:
        splog.debug('psutil is not installed... not reporting open files')
        return set()
    process = psutil.Process(os.getpid())
    temp = set([f.path for f in process.open_files()])
    if prv is not None:
        leaked_files = temp - prv
        if leaked_files:
            splog.debug("Leaked files:"+', '.join( leaked_files))
        else: splog.debug("No file leaks detected.")
        splog.debug("Open files"+', '.join(prv))
    return temp

def start():
    tracemalloc.start()
    return

def stop():
    tracemalloc.stop()
    return
    
def clean(prv=None):
    pltclose('all')
    gc.collect()
    try:
        for obj in gc.get_objects():
            if isinstance(obj, fitsHDUList): obj.close()
        for obj in gc.get_objects():
            if isinstance(obj, io.IOBase) and not obj.closed:
                if hasattr(obj, "name") and isinstance(obj.name, str) and (obj.name.endswith(".fits") or obj.name.endswith(".fits.gz")):
                    splog.debug(f"Closing: {obj.name}")
                    obj.close()
    except Exception as e:
        splog.debug(f'{e}')
    try:
        all_objects = gc.get_objects()
        # Filter out the ones that are instances of the fitsrec class
        fitsrec_objects = [obj for obj in all_objects if (isinstance(obj, FITS_record) or isinstance(obj,FITS_rec))]
        if len(fitsrec_objects) > 0:
            splog.debug(','.join(fitsrec_objects))
        for obj in fitsrec_objects:
            del obj
        gc.collect()
    except Exception as e:
        splog.debug(f"{e}")

    try:
        after_gc = check(prv=prv)
    except:
        pass
    return
