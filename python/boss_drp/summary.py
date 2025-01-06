from boss_drp import idlspec2d_dir
import os.path as ptt
from boss_drp.field import Field, field_to_string
from os import makedirs, getenv
from typing import Optional

class Summary_names:
    def __init__(self):#, #indir: str, run2d: str, outroot: Optional[str] = None,
                # field: Optional[str] = None, mjd: Optional[int] = None,
               #  dev: bool = False, epoch: bool = False, custom: Optional[str] = None,
               #  allsky: bool = False, tmpext: str = '', outdir: Optional[str] = None):
        
        # Initialize attributes to None to prevent AttributeError
        self.outdir = None
        self.spAllfile = None
        self.spAlllitefile = None
        self.splinefile = None
        self.spAlldatfile = None
        self.datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spall_dm.par')
        self.line_datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spzline_dm.par')
        
    def set(self, indir: str, run2d: str, outroot: Optional[str] = None,
                 field: Optional[str] = None, mjd: Optional[int] = None,
                 dev: bool = False, epoch: bool = False, custom: Optional[str] = None,
                 allsky: bool = False, tmpext: str = '', outdir: Optional[str] = None):
                 
        self.build(indir, run2d, outroot=outroot, field=field, mjd=mjd, dev=dev,
                    epoch=epoch, custom=custom, allsky=allsky, outdir=outdir)



        # Initialize the BK (backup) nested object
        self.bk = self.BK(self)
        self.temp = self.TEMP(self,tmpext)

    class BK:
        def __init__(self, fnames):
            # Create modified versions of the parent Fnames filenames
            self.fnames = fnames
            self.set()

        def modify_filename(self, filename):
            # Define the modification logic (e.g., add "_bk" before the file extension)
            return f"{filename.replace('spAll','bkup/spAll')}.bkup"
        
        def set(self, flag=None):
            self.spAllfile = self.modify_filename(self.fnames.spAllfile)+f'-{flag}'
            self.spAlllitefile = self.modify_filename(self.fnames.spAlllitefile)+f'-{flag}'
            self.splinefile = self.modify_filename(self.fnames.splinefile)+f'-{flag}'
            self.spAlldatfile = self.modify_filename(self.fnames.spAlldatfile)+f'-{flag}'
            makedirs(ptt.dirname(self.spAllfile), exist_ok = True)

    class TEMP:
        def __init__(self, fnames, tmpext):
            # Create modified versions of the parent Fnames filenames
            self.fnames = fnames
            self.tmpext = tmpext
            self.spAllfile = self.modify_filename(fnames.spAllfile, self.tmpext)
            self.spAlllitefile = self.modify_filename(fnames.spAlllitefile, self.tmpext)
            self.splinefile = self.modify_filename(fnames.splinefile, self.tmpext)
            self.spAlldatfile = self.modify_filename(fnames.spAlldatfile, self.tmpext)

        @staticmethod
        def modify_filename(filename, tmpext):
            # Define the modification logic (e.g., add "_bk" before the file extension)
            return f"{filename.replace('.gz',tmpext+'.gz')}"
    

    def build(self, indir, run2d, outroot=None, field=None, mjd=None, dev=False,
                    epoch=False, custom=None, allsky=False, outdir=None):
                
        if outroot is not None:
            self.spAllfile     = ptt.join(outroot+'.fits.gz')
            self.spAlllitefile = ptt.join(outroot+'-lite'+'.fits.gz')
            self.splinefile    = ptt.join(outroot+'Line'+'.fits.gz')
            self.spAlldatfile  = ptt.join(outroot+'.dat.gz')
        else:
            cc = False
            if custom is not None:
                cc= True
            elif field is not None:
                field = field_to_string(field)
            if field is not None and mjd is not None:
                field_class = Field(indir, run2d, field, custom_name=custom)
                specfull_dir  = field_class.spec_dir(mjd,epoch=epoch)
                mjd = str(mjd)
                self.spAllfile     = ptt.join(specfull_dir, 'spAll-'+field+'-'+mjd+'.fits.gz')
                self.spAlllitefile = ptt.join(specfull_dir, 'spAll-lite-'+field+'-'+mjd+'.fits.gz')
                self.splinefile    = ptt.join(specfull_dir, 'spAllLine-'+field+'-'+mjd+'.fits.gz')
                self.spAlldatfile  = ptt.join(specfull_dir, 'spAll-'+field+'-'+mjd+'.dat.gz')
            else:
                fflags = []
                if run2d is not None:
                    spall_dir = Summary_dir(indir, run2d, epoch=epoch, custom_name=custom)
                else:
                    spall_dir = Summary_dir(indir, '', epoch=epoch, custom_name=custom)
                if outdir is not None:
                    spall_dir = outdir

                fflags = [run2d if run2d else '', custom if custom else '', 'epoch' if epoch else '']
                fflags = '-'.join(filter(None, fflags))  # Remove empty strings
                if len(fflags) > 0:
                    fflags = '-'+fflags
                self.spAllfile     = ptt.join(spall_dir, 'spAll'+fflags+'.fits.gz')
                self.spAlllitefile = ptt.join(spall_dir, 'spAll-lite'+fflags+'.fits.gz')
                self.splinefile    = ptt.join(spall_dir, 'spAllLine'+fflags+'.fits.gz')
                self.spAlldatfile  = ptt.join(spall_dir, 'spAll'+fflags+'.dat.gz')

        if dev:
            self.spAllfile = self.spAllfile.replace('spAll','spAll_dev')
            self.spAlllitefile = self.spAlllitefile.replace('spAll','spAll_dev')
            self.splinefile = self.splinefile.replace('spAllLine','spAllLine_dev')
            self.spAlldatfile = self.spAlldatfile.replace('spAll','spAll_dev')
        self.outdir = ptt.dirname(self.spAllfile)
        return

summary_names = Summary_names()

class FieldList_name:
    def __init__(self):
        self.outdir = None
        self.run2d = None
        self.epoch = False
        self.topdir = None
        self.custom = False
        self.custom_name = None
        self.tmpext = '.tmp'
        self.html = {}
        self.logfile = None
        self.basehtml = None
    def build(self, topdir: str, run2d: str, epoch: bool = None,
                 outdir: Optional[str] = None,
                 custom_name: Optional[str] = None,
                 logfile: Optional[str] = None,
                 tmpext: Optional[str] = None):
                 
        if outdir is not None:
            self.outdir = outdir
        if run2d is not None:
            self.run2d = run2d
        if epoch is not None:
            self.epoch = epoch
        if topdir is not None:
            self.topdir = topdir
        if custom_name is not None:
            custom = True
            self.custom = custom
            self.custom_name = custom_name
        if tmpext is not None:
            self.tmpext = tmpext
            
        if self.outdir is None:
            self.outdir = Summary_dir(self.topdir, self.run2d, epoch=False, outdir=None, custom=False, custom_name=None)
        if self.epoch:
            self.name = ptt.join(self.outdir, 'fieldlist-'+self.run2d+'-epoch.fits')
        elif self.custom:
            self.name = ptt.join(self.outdir, 'fieldlist-'+self.run2d+'-'+self.custom_name+'.fits')
        else:
            self.name = ptt.join(self.outdir, 'fieldlist-'+self.run2d+'.fits')

        self.html = {}
        self.html['fieldlist'] = 'fieldlist{obs}.html'
        self.html['fieldlist_mjd'] = 'fieldlist{obs}-mjdsort.html'
        self.html['fieldquality'] = 'fieldquality{obs}.html'
        self.html['fieldquality_mjd'] =  'fieldquality{obs}-mjdsort.html'

        if logfile is None:
            self.tmpext = '.tmp'
            self.logfile = ptt.join(self.outdir, self.name.replace('.fits','.log'))
        else:
            self.tmpext = ptt.basename(logfile).replace('.log','.tmp').replace('fieldlist-','')
            self.logfile = logfile

fieldlist_name = FieldList_name()

def Base_dir(topdir, run2d):
    return ptt.join(topdir, run2d)

def Summary_dir(topdir, run2d, epoch=False, custom_name=None, custom=False, outdir=None):
    if custom_name is not None:
        custom = True
    if outdir is None:
        if epoch:
            outdir = ptt.join(topdir, run2d, 'summary','epoch')
        elif custom:
            outdir = ptt.join(topdir, run2d, 'summary',custom_name)
        else:
            outdir = ptt.join(topdir, run2d, 'summary','daily')
    makedirs(outdir, exist_ok = True)
    return outdir
