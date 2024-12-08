from os import getenv

def chpc2html(fpath):
    sas_base_dir = getenv('SAS_BASE_DIR', default=None)
    if sas_base_dir is None:
        sas_base_dir = "/uufs/chpc.utah.edu/common/home/sdss50"
    return(fpath.replace(sas_base_dir,"https://data.sdss5.org/sas/"))
