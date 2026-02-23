
from os import getenv
import re
from sdss_access import Path
import os
path = Path(release = os.getenv('TREE_VER',None), preserve_envvars=True) 

def chpc2html(fpath):
    sas_base_dir = getenv('SAS_BASE_DIR')
    if sas_base_dir is None:
        # NOTE: no number here — just sdss
        sas_base_dir = "/uufs/chpc.utah.edu/common/home/sdss"
        
    # Build a regex that matches sas_base_dir plus digits
    pattern = re.escape(sas_base_dir) + r"\d*(?=/|$)"

    
    # Replacement URL base
    base_url = path.remote_base.rstrip('/') + '/sas/'
    
    # Perform replacement
    url = re.sub(pattern, base_url, fpath)
    
    # Split protocol
    if url.startswith("https://"):
        proto = "https://"
        rest = url[len("https://"):]
    elif url.startswith("http://"):
        proto = "http://"
        rest = url[len("http://"):]
    else:
        proto = ""
        rest = url
        
    # Collapse accidental repeated slashes AFTER protocol
    rest = re.sub(r'/+', '/', rest)
    
    return proto + rest
