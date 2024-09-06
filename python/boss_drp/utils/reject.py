
from astropy.io import fits
import os.path as ptt

class Reject:
    def __init__(self, frame, hdr):
        self.frame = frame
        self.hdr = hdr
        self.ffs = 0
        self.ff = 0
        self.ne = 0
        self.hgcd = 0
        self.hear = 0
        self.obs  = ''
        self.favor= ''
        self.hartmann = ''
        self._lamps()
        
    def _lamps(self):
        self.ffs  = self.hdr.get('FFS', '0 0 0 0 0 0 0 0').split().count('1')
        self.ff   = self.hdr.get('FF',  '0 0 0 0').split().count('1')
        self.ne   = self.hdr.get('NE',  '0 0 0 0').split().count('1')
        self.hgcd = self.hdr.get('HGCD','0 0 0 0').split().count('1')
        self.hear = self.hdr.get('HEAR','0 0 0 0').split().count('1')
        self.flavor = self.hdr.get('FLAVOR', '')
        cart = str(self.hdr.get('CARTID', 'FPS-N'))
        self.hartmann = self.hdr.get('HARTMANN', 'out').lower()
        if 'FPS-N' in cart:
            self.obs = 'APO'
        elif 'FPS-S' in cart:
            self.obs = 'LCO'
        else:
            self.obs = 'APO'
    

    def check(self, splog, hartmann = False):
        if self.flavor.lower() in ['arc','calibration']:
            return self.check_arc(splog, hartmann=hartmann)
        elif self.flavor.lower() == 'flat':
            return self.check_flat(splog)
        elif self.flavor.lower() in ['target', 'science']:
            return self.check_science(splog)
        splog.info('Frame is not Arc, Flat, or Science')
        return True


    def check_flat(self, splog):
        if self.obs == 'APO':
            if self.ne > 0:
                splog.info(f'Warning: Reject Flat: {self.ne}/{4} Ne lamps are On! ({ptt.basename(self.frame)})')
                return True
            elif self.hgcd > 0:
                splog.info(f'Warning: Reject Flat: {self.hgcd}/{4} HgCd lamps are On! ({ptt.basename(self.frame)})')
                return True
            else:
                pass
        elif self.obs == 'LCO':
            if self.ne > 0:
                splog.info(f'Warning: Reject Flat: {self.ne}/{4} Ne lamps are On! ({ptt.basename(self.frame)})')
                return True
            elif self.hgcd > 0:
                splog.info(f'Warning: Reject Flat: self.{hgcd}/{4} HeAr lamps are On! ({ptt.basename(self.frame)})')
                return True
            else:
                pass
        if self.ff < 4:
            splog.info(f'Warning: {4-self.ff}/4 Flat-field lamps turned off ({ptt.basename(self.frame)})')

        if self.ffs < 8:
            splog.info(f'Warning: Reject Flat: {self.ffs}/8 Flat-field screens not closed! ({ptt.basename(self.frame)})')
            return True

        if 'out' not in self.hartmann:
            splog.info(f'Warning Hartmann doors closed ({ptt.basename(self.frame)})')
            return True
        
        return False

    def check_arc(self, splog, hartmann=False):
        if self.obs == 'APO':
            if self.ne < 4 and self.hgcd < 4:
                splog.info(f'WARNING: Reject arc: Neither Ne nor HgCd lamps are On! ({ptt.basename(self.frame)})')
                return True
            elif self.ne < 4:
                splog.info(f'Warning: {4-self.ne}/{4} Ne lamps are off ({ptt.basename(self.frame)})')
            elif self.hgcd < 4:
                splog.info(f'Warning: {4-self.hgcd}/{4} HgCd lamps are off ({ptt.basename(self.frame)})')
            else:
                pass
        elif self.obs == 'LCO':
            if self.ne < 4 and self.hear < 4:
                splog.info('WARNING: Reject arc: Neither Ne nor HeAr lamps are off! ({ptt.basename(self.frame)})')
                return True
            elif self.ne < 4:
                splog.info(f'Warning: {4-self.ne}/{4} Ne lamps are off ({ptt.basename(self.frame)})')
            elif self.hear < 4:
                splog.info(f'Warning: {4-self.hear}/{4} HeAr lamps are off ({ptt.basename(self.frame)})')
            else:
                pass
        if self.ff > 0:
            splog.info(f'Warning: Reject Arc: {self.ff}/4 Flat-field lamps turned on! ({ptt.basename(self.frame)})')
            return True
      
        if self.ffs < 8:
            splog.info(f'Warning: Reject Arc: {self.ffs}/8 Flat-field screens not closed! ({ptt.basename(self.frame)})')
            return True

        if not hartmann:
            if 'out' not in self.hartmann:
                splog.info(f'Warning Hartmann doors closed ({ptt.basename(self.frame)})')
                return True
        
        return False

    def check_science(self, splog):
        if self.obs == 'APO':
            if self.ne > 0:
                splog.info(f'Warning: Reject Science: {self.ne}/{4} Ne lamps are On! ({ptt.basename(self.frame)})')
                return True
            elif self.hgcd > 0:
                splog.info(f'Warning: Reject Science: {self.hgcd}/{4} HgCd lamps are On! ({ptt.basename(self.frame)})')
                return True
            else:
                pass
        elif self.obs == 'LCO':
            if self.ne > 0:
                splog.info(f'Warning: Reject Science: {self.ne}/{4} Ne lamps are On! ({ptt.basename(self.frame)})')
                return True
            elif self.hgcd > 0:
                splog.info(f'Warning: Reject Science: self.{hgcd}/{4} HeAr lamps are On! ({ptt.basename(self.frame)})')
                return True
            else:
                pass
        if self.ff > 0:
            splog.info(f'Warning: Reject Science: {self.ff}/4 Flat-field lamps turned on! ({ptt.basename(self.frame)})')
            return True
        if self.ffs > 0:
            splog.info(f'Warning: Reject Science: {self.ffs}/8 Flat-field screens closed! ({ptt.basename(self.frame)})')
            return True
        return False
