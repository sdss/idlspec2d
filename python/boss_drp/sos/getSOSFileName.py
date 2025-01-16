from boss_drp.utils import getcard

import os.path as ptt
from astropy.io import fits

def remove_all_extensions(filename):
    while True:
        base, ext = ptt.splitext(filename)
        if ext:  # If there's an extension, keep removing
            filename = base
        else:    # Stop when no extension is left
            break
    return filename


def getSOSFileName(filename):
    hdr = fits.header(filename)
    confID = getcard(hdr, 'CONFID', default=0, noNaN=True)
    confID = f"{confID:06}"
    fieldid = getcard(hdr, 'FIELDID', default=0, noNaN=True)

    
    filename = ptt.basename(filename)
    filename = remove_all_extensions(filename)
    _, camera, exposure = filename.split('-') #root (sdR), camera name, exposure number
    
    return f"sci-{confID:06}-{camera}-{exposure}.fits"
