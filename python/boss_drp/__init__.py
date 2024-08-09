#import subprocess
#import os
#from pkg_resources import resource_filename

from sdsstools import get_package_version
#from sdsstools import get_config, Configuration
import os
import numpy as np

__version__ = get_package_version(__file__, 'boss_drp') or 'dev'

try:
    daily_dir = os.getenv('DAILY_DIR')
except:
    daily_dir = None
    pass
if daily_dir is None:
    daily_dir = os.path.join(os.getenv('HOME'),'daily')
    os.environ['DAILY_DIR'] = daily_dir
    
