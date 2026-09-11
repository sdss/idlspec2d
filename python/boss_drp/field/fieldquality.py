from boss_drp import idlspec2d_dir
from boss_drp.field.generations import generations
from pathlib import Path
import numpy as np
from ruamel.yaml import YAML
yaml = YAML()
### Field Quality Limits

class modeFieldQuality:
    def __init__(self, blue, red, exptime = 3600, threshold_scale=1, design_completion=False, mag15=False):
        """
        Defines the FIELDQUALITY Limits per designmode
            paramters:
            ----------
            blue: float
                blue Field SN2 threshold
            red:  float
                red Field SN2 threshold
            exptime: float
                the exposure time to scale the measure SN2 to match
            threshold_scale: float
                Scale Factor to scale the Field SN2 requirements to the full pipeline values
            design_completion: bool
                Flag indicating if design completion criteria rather then epoch criteria
            mag15: bool
                Flag indicating if criteria uses Mag15 SN2
        """
        
        self.blue = blue*threshold_scale
        self.red = red*threshold_scale
        self.exptime = exptime
        self.design_completion = design_completion
        self.mag15 = mag15

    def check(self, min_sn2_b, min_sn2_r, exptime =  None):
        """
        Checks the SN2 values for the given criteria
        
        Parameters:
        -----------
        min_sn2_b: float
            the SN2 blue value (minimum of the b1 and b2)
        min_sn2_r: float
            the SN2 red value (minimum of the r1 and r2)
        exptime: float
            the observed exposure time for scaling (if set)

        Returns:
        --------
        iqual: int
            the incoming quality flag (associated with ['bad', 'marginal', 'good'])
        """
        if exptime:
            min_sn2_b=min_sn2_b/(exptime/self.exptime)
            min_sn2_r=min_sn2_r/(exptime/self.exptime)

        if np.isnan(self.red):
            if np.isnan(min_sn2_b):
                return 0
            if (min_sn2_b < self.blue):
                return 0
        if np.isnan(self.blue):
            if np.isnan(min_sn2_r):
                return 0
            if (min_sn2_r < self.red):
                return 0
        if np.isnan(min_sn2_r) and np.isnan(min_sn2_b):
            return 0
        if ((min_sn2_b < self.blue) or (min_sn2_r < self.red)):
            return 0

        return 2

class FieldQuality:
    def __init__(self):
        with open(Path(idlspec2d_dir) / 'python'/'boss_drp'/'etc' / 'fieldquality.yml') as f:
            requirements = yaml.load(f)
        requirements = yaml.load(Path(idlspec2d_dir) / 'python'/'boss_drp'/'etc' / 'fieldquality.yml')
        mode_to_qual = requirements['mode_to_qual']
        thresholds = requirements['Qual_thresholds']
        _defaults = {'blue': np.nan, 'red': np.nan, 'exptime': 3600.0, 
                     'threshold_scale': 1.0, 'design_completion': False, 
                     'mag15': False}
        self.modes = {}
        for mode, limit_type in mode_to_qual.items():
            limits = thresholds[limit_type]
            for key in _defaults:
                if key not in limits:
                    limits[key] = _defaults[key]
                elif limits[key] is None:
                    limits[key] = _defaults[key]
            self.modes[mode.lower()] = modeFieldQuality(limits['blue'], limits['red'], 
                                                        exptime=limits['exptime'], 
                                                        threshold_scale=limits['threshold_scale'],
                                                        design_completion=limits['design_completion'])

    
    def check(self, row, min_sn2_b, min_sn2_r, min_sn2_15_b=None, min_sn2_15_r=None, daily = False):
        """
        Determine the FieldQuality Criteria and check it
        
        Parameters:
        -----------
        row: Table.row
            astropy Table row of the fieldlist
        min_sn2_b: float
            the SN2 blue value (minimum of the b1 and b2)
        min_sn2_r: float
            the SN2 red value (minimum of the r1 and r2)
        min_sn2_15_b: float
            the SN2 blue value (minimum of the b1 and b2) at Mag = 15 
        min_sn2_15_r: float
            the SN2 red value (minimum of the r1 and r2) at Mag = 15 
        daily: Bool
            set if daily coadds
            
        return:
        -------
        iqual: int
            the incoming quality flag (associated with ['bad', 'marginal', 'good'])

        """
        
        #if we want to scale the criteria to number of exposure... however, that could be misleading
        #scale = True if daily else False
        scale = False
        exptime = row['EXPTIME'] if scale else None # to scale to number of exposures
        
        mode = row['DESIGN_MODE'].lower()
        prog = row['PROGRAMNAME'].strip() if 'PROGRAMNAME' in row else ''

        row_gen = generations.find_generation(row['FIELD'], row['MJD'], row['OBSERVATORY'])
        
        if (row_gen is None):
            # FieldID/MJD is not in any generation, so we cannot determine the quality
            return 0
        row_gen = [r.lower() for r in row_gen]
        if 'legacy' in row_gen:
            # Pre SDSS-V Plates
            mode = 'legacy' if int(row['MJD']) > 58029 else 'legacy_early'
            if prog.upper() in ['ELG_NGC','ELG_SGC']:
                #--- JEB 2018-05-23: if elg plate, plate is 'good' no matter what SN2
                return 2
        elif 'plates' in row_gen:
            # SDSS-V Plates
            mode = 'plates'
            if 'RM' in prog:
                mode = 'plates_rm'
            mode = f'sdss5_{mode}'
        elif ('fps' in row_gen) and ('sdssv' in row_gen):
            # SDSS-V FPS
            mode = f'sdss5_{mode}'
            mode = mode.replace('_no_apogee_skies','')
        elif ('fps' in row_gen) and ('as5' in row_gen):
            # AS5 FPS
            mode = f'as5_{mode}'

        if self.modes[mode].design_completion:
            """ dark_plane only has a design completion criteria so determine number
                of designs and use that as the scale factor """
            exptime =  len(row['DESIGNS'].split(' ')) * 900.0

        if self.modes[mode].mag15:
            min_sn2_b = min_sn2_15_b
            min_sn2_r = min_sn2_15_r

        iqual = self.modes[mode].check(min_sn2_b, min_sn2_r, exptime = exptime)
        return min(2, iqual)


fieldquality = FieldQuality()
