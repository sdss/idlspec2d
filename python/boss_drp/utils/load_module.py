import os.path as ptt
from os import getenv
import importlib
import sys


def import_module(path = '/uufs/chpc.utah.edu/sys/pkg/modules/init/python3'):
    spec = importlib.util.spec_from_loader(ptt.basename(path),
                importlib.machinery.SourceFileLoader(ptt.basename(path), path))
    pymodule = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pymodule)
    sys.modules[ptt.basename(path)] = pymodule
    return(pymodule)

def load_env(key, default=None):
    val = getenv(key)
    if val is None:
        val = default
        if val is None:
            print('ERROR: '+key+' is not set')
            exit()
    return(val)


class ModuleException(Exception):
    pass


def load_module():
    try:
        mod = __import_module(path = '/uufs/chpc.utah.edu/sys/pkg/modules/init/python3')
    except Exception as e:
        raise ModuleException('Can not access the Utah CHPC modules') from None
    return(mod.module)
