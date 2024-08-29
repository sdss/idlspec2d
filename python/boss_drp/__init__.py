#import subprocess
#import os
#from pkg_resources import resource_filename

from sdsstools import get_package_version
#from sdsstools import get_config, Configuration
import os
import numpy as np
import warnings
__version__ = get_package_version(__file__, 'boss_drp') or 'dev'




try:
    daily_dir = os.getenv('DAILY_DIR')
except:
    daily_dir = None
    pass
if daily_dir is None:
    daily_dir = os.path.join(os.getenv('HOME'),'daily')
    os.environ['DAILY_DIR'] = daily_dir

class MissingEnvVarWarning(UserWarning):
    pass

try:
    idlspec2d_dir = os.getenv('IDLSPEC2D_DIR')
except:
    idlspec2d_dir = None
if idlspec2d_dir is None:
    warnings.warn('IDLSPEC2D_DIR ENV Variable is not set',MissingEnvVarWarning)
    os.environ['IDLSPEC2D_DIR']  = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
idlspec2d_dir = os.getenv('IDLSPEC2D_DIR')



favicon ="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png"
