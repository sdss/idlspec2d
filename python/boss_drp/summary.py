from boss_drp import idlspec2d_dir
import os.path as ptt
from boss_drp.field import Field, field_to_string
from os import makedirs, getenv
from typing import Optional
from pathlib import Path

def _f2pq(fpath):
        fpath = Path(fpath)
        return str(Path(str(fpath).removesuffix(''.join(fpath.suffixes))).with_suffix('.parquet'))

class Summary_names:
    def __init__(self):
        
        # Initialize attributes to None to prevent AttributeError
        self.outdir = None
        self.spAllfile = None
        self.spAlllitefile = None
        self.splinefile = None
        self.spcalibfile = None
        self.epoch = False
        self.custom = None
        self.allsky = False
        self.fmjd_ver = False
        self.MJD_dir = None
        self.spAllfile_parquet=None
        self.spAlllitefile_parquet = None
        self.splinefile_parquet = None
        self.spcalibfile_parquet = None
        self.daily_spAll_parquet = None
        self.daily_spline_parquet = None
        self._build_args = {}

        self.datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spall_dm.par')
        self.line_datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spzline_dm.par')
        self.spcalib_datamodel = ptt.join(idlspec2d_dir, 'datamodel', 'spcalib_dm.par')
        
    def set(self, indir: str, run2d: str, outroot: Optional[str] = None,
                 field: Optional[str] = None, mjd: Optional[int] = None,
                 dev: bool = False, epoch: bool = False, custom: Optional[str] = None,
                 allsky: bool = False, tmpext: str = '', outdir: Optional[str] = None,
                 MJD_dir: Optional[str] = None, obs: Optional[str] = None):

        self._build_args = locals()  # Store the arguments for potential future use
        self._build_args.pop('self')  # Remove 'self' from the stored arguments

        self.build(indir, run2d, outroot=outroot, field=field, mjd=mjd, dev=dev, MJD_dir=MJD_dir,
                    epoch=epoch, custom=custom, allsky=allsky, outdir=outdir, obs=obs)

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
            self.spcalibfile = self.modify_filename(self.fnames.spcalibfile)+f'-{flag}'

            self.spAllfile_parquet = self.modify_filename(self.fnames.spAllfile_parquet)+f'-{flag}'
            self.spAlllitefile_parquet = self.modify_filename(self.fnames.spAlllitefile_parquet)+f'-{flag}'
            self.splinefile_parquet = self.modify_filename(self.fnames.splinefile_parquet)+f'-{flag}'
            self.spcalibfile_parquet = self.modify_filename(self.fnames.spcalibfile_parquet)+f'-{flag}'

        def mkdir(self):
            makedirs(ptt.dirname(self.spAllfile), exist_ok = True)

    class TEMP:
        def __init__(self, fnames, tmpext):
            # Create modified versions of the parent Fnames filenames
            self.fnames = fnames
            self.tmpext = tmpext
            self.spAllfile = self.modify_filename(fnames.spAllfile, self.tmpext)
            self.spAlllitefile = self.modify_filename(fnames.spAlllitefile, self.tmpext)
            self.splinefile = self.modify_filename(fnames.splinefile, self.tmpext)
            self.spcalibfile  = self.modify_filename(fnames.spcalibfile, self.tmpext)

            self.spAllfile_parquet = self.modify_filename(fnames.spAllfile_parquet, self.tmpext)
            self.spAlllitefile_parquet = self.modify_filename(fnames.spAlllitefile_parquet, self.tmpext)
            self.splinefile_parquet = self.modify_filename(fnames.splinefile_parquet, self.tmpext)
            self.spcalibfile_parquet  = self.modify_filename(fnames.spcalibfile_parquet, self.tmpext)

        @staticmethod
        def modify_filename(filename, tmpext):
            # Define the modification logic (e.g., add "_bk" before the file extension)
            
            if tmpext == '': return filename
            p = Path(filename)
            p = p.with_name(f"{tmpext}_{p.name}")
            return str(p)

            # return f"{filename.replace('.gz',tmpext+'.gz')}"
    

    def build(self, indir: str, run2d: str, outroot: Optional[str] = None,
              field: Optional[str] = None, mjd: Optional[int] = None,
              dev: bool = False, epoch: bool = False, custom: Optional[str] = None,
              allsky: bool = False,  outdir: Optional[str] = None,
              MJD_dir: Optional[str] = None, obs: Optional[str] = None):
                        
        self.epoch = epoch
        self.custom = custom
        self.allsky = allsky
        self.fmjd_ver = False
        self.MJD_dir = MJD_dir
        if len(self._build_args) == 0:
            self._build_args = locals()  # Store the arguments for potential future use
            self._build_args.pop('self')  # Remove 'self' from the stored arguments

        if outroot is not None:
            self.spAllfile     = ptt.join(outroot+'.fits.gz')
            self.spAlllitefile = ptt.join(outroot+'-lite'+'.fits.gz')
            self.splinefile    = ptt.join(outroot+'Line'+'.fits.gz')
            self.spcalibfile   = ptt.join(outroot+'-calib_qa.fits')
        else:
            cc = False
            if custom is not None:
                cc= True
            elif field is not None:
                field = field_to_string(field)
            if field is not None and mjd is not None:
                field_class = Field(indir, run2d, field, custom_name=custom)
                spall_dir  = field_class.spec_dir(mjd,epoch=epoch)
                mjd = str(mjd)
                self.fmjd_ver = True
                fflags = [f"{field}-{mjd}"]
            else:
                fflags = [f'{run2d}']
                if run2d is not None:
                    spall_dir = Summary_dir(indir, run2d, epoch=epoch, custom_name=custom)
                else:
                    spall_dir = Summary_dir(indir, '', epoch=epoch, custom_name=custom)
                if outdir is not None:
                    spall_dir = outdir
                if cc:
                    fflags.append(custom)
                if epoch:
                    fflags.append('epoch')

            fflags = f'-{"-".join(fflags)}' if len(fflags) > 0 else ''
            self.spAllfile     = ptt.join(spall_dir, 'spAll'+fflags+'.fits.gz')
            self.spAlllitefile = ptt.join(spall_dir, 'spAll-lite'+fflags+'.fits.gz')
            self.splinefile    = ptt.join(spall_dir, 'spAllLine'+fflags+'.fits.gz')
            self.spcalibfile   = ptt.join(spall_dir, 'spCalib_QA'+fflags+'.fits')

        if dev:
            self.spAllfile = self.spAllfile.replace('spAll','spAll_dev')
            self.spAlllitefile = self.spAlllitefile.replace('spAll','spAll_dev')
            self.splinefile = self.splinefile.replace('spAllLine','spAllLine_dev')
            self.spcalibfile  = self.spcalibfile.replace('spCalib_QA','spCalib_QA_dev')

        self.spAllfile_parquet     = _f2pq(self.spAllfile)
        self.spAlllitefile_parquet = _f2pq(self.spAlllitefile)
        self.splinefile_parquet    = _f2pq(self.splinefile)
        self.spcalibfile_parquet   = _f2pq(self.spcalibfile)

        self.daily_spAll_parquet   = 'spAll-{run2d}-{mjd}_{obs}.parquet'
        self.daily_spline_parquet   = 'spAllLine-{run2d}-{mjd}_{obs}.parquet'

        self.outdir = ptt.dirname(self.spAllfile)

        if self.MJD_dir is None:
            self.MJD_dir = getenv('BOSS_SPECTRO_SCRATCH', None)
        if (self.MJD_dir is None) or (self.MJD_dir == ''):
            self.MJD_dir = outdir
        else:
            redux = getenv('BOSS_SPECTRO_REDUX')
            scratch_dir = Path(self.MJD_dir)
            try:
                self.MJD_dir = str(scratch_dir / Path(self.outdir).relative_to(redux))
            except:
                pass
        self.MJD_dir = Path(self.MJD_dir) / 'mjd'
        self.MJD_dir = str(self.MJD_dir)

        return


    def clone(self, **overrides):
        args = self._build_args.copy()
        args.update(overrides)

        new = Summary_names()
        new.set(**args)
        new.bk = new.BK(new)
        new.temp = new.TEMP(new, "")
        return new

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
            self.outdir = Summary_dir(self.topdir, self.run2d, epoch=self.epoch, outdir=None, custom=self.custom, custom_name=self.custom_name)
        if self.epoch:
            self.name = ptt.join(self.outdir, 'fieldlist-'+self.run2d+'-epoch.fits')
        elif self.custom:
            self.name = ptt.join(self.outdir, 'fieldlist-'+self.run2d+'-'+self.custom_name+'.fits')
        else:
            self.name = ptt.join(self.outdir, 'fieldlist-'+self.run2d+'.fits')
        self.parquet = _f2pq(self.name)

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
        
    def temp(self,parquet=False):
        if not parquet:
            return self.name.replace('fieldlist-', 'tmp_fieldlist-')
        return self.parquet.replace('fieldlist-', 'tmp_fieldlist-')



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
