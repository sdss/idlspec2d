from .field_to_string import field_to_string
from .fieldgroup import fieldgroup
from .Fieldtype import Fieldtype
#!/usr/bin/env python3
import numpy as np
import os.path as ptt
from boss_drp.utils.splog import splog

class Field:
    def __init__(self, topdir, run2d, field, mjd=None,custom=False, custom_name=None,
                epoch = False, run1d = None, obs = None):
        self.topdir = topdir
        self.run2d  = run2d
        self.run1d  = run1d
        self.obs    = obs
        self.mjd    = mjd
        self.field  = field
        self.plan2d = None
        self.plancomb = None
        self.spField = None
        self.configIDs = None
        self.custom = custom
        self._custom = custom
        if custom_name is not None:
            self.custom = True
        self.custom_name = custom_name
        self.epoch = epoch
        self.set()
        
    def set(self):
        if not self.custom:
            self.field_str = field_to_string(self.field)
        else:
            self.field_str = self.field
        if self.custom_name is not None:
            self.custom = True
        else:
            self.custom = self._custom
        
        self.fieldgroup = self.setgroup()
        try:
            self.type = Fieldtype(fieldid = self.field, mjd = self.mjd)
        except:
            if self.custom is False and self.field != '*':
                splog.warning('Undetermined Field.type')
            self.type = Fieldtype()
            
    def setgroup(self):
        return fieldgroup(self.field, custom=self.custom)
#        field = str(self.field)
#        if self.custom:
#            if len(field.split('_')) == 1:
#                return(field)
#            else:
#                return('_'.join(field.split('_')[:-1]))
#        elif field.isnumeric():
#            return(str(np.floor(int(field)/1000).astype(int)).zfill(3)+'XXX')
#        elif field == '*':
#            zfield = field_to_string(0)
#            return(fieldgroup(zfield).replace('0','?'))
#       else:
#            return(field)
    
    def spec1d_dir(self, run1d=None, epoch =None):
        if run1d is None:
            run1d = self.run1d
        return ptt.join(self.dir(epoch = epoch), run1d)

    def dir(self, epoch=None):
        if epoch is None: epoch = self.epoch
        topdir2d = ptt.join(self.topdir, self.run2d)
        if self.custom:
            dir_ = ptt.join(topdir2d, 'fields', self.fieldgroup, self.custom_name)
        else:
            dir_ = ptt.join(topdir2d, 'fields', self.fieldgroup, self.field)
        if epoch is True:
            dir_ = ptt.join(dir_, 'epoch')
        return(dir_)
    
    
    def spec_dir(self, mjd=None, epoch=None, full = True, pathbase=None):
        if epoch is None: epoch = self.epoch
        if mjd is None: mjd = self.mjd
        mjd = str(mjd)
        stype = 'full' if full else 'lite'
        fs = self.field_str
        if epoch: dirtype = 'epoch'
        elif self.custom_name is not None: dirtype = self.custom_name
        else: dirtype = 'daily'
        if pathbase is None:
            pathbase = ptt.join(self.topdir)
        _base = ptt.join(pathbase, self.run2d, 'spectra')
        
        return ptt.join(_base, dirtype, stype, self.fieldgroup, fs, mjd)
        
    def png_dir(self, run1d=None, mjd=None, epoch=None, pathbase=None):
        if run1d is None:
            run1d = self.run1d
        if mjd is None:
            mjd = self.mjd
        if pathbase is None:
            pathbase = ptt.join(self.topdir)
        if epoch is None: epoch = self.epoch
        mjd = str(mjd)
        if epoch: dirtype = 'epoch'
        elif self.custom_name is not None: dirtype = self.custom_name
        else: dirtype = 'daily'
        
        imagebase = ptt.join(pathbase, self.run2d, 'images', dirtype, run1d)

        return ptt.join(imagebase, self.fieldgroup, self.field_str, self.field_str+'-'+mjd)
