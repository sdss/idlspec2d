import boss_drp
import json
from pathlib import Path

class Generations:
    def __init__(self):
        self.mapping = None
        self._map_file = Path(boss_drp.__file__).parent / 'etc' / 'generations.json'
    def load(self):
        """
        Load the generations mapping from the JSON file 
        """
        if self.mapping is None:
            with open(self._map_file , "r") as f:
                self.mapping = json.load(f)

    def get(self, gen, field = False, mjd = False, obs= None):
        """
        Returns the field or mjd range for the given generation and observatory
        Parameters:
        -----------
        gen: str
            the generation to get the range for
        field: bool
            if True, return the field range
        mjd: bool
            if True, return the mjd range
        obs: str
            the observatory to get the range for (apo, lco)
        Returns:
        -------- 
        range: list
            the field or mjd range for the given generation and observatory
        """
        if self.mapping is None:
            self.load()
        if field:
            return self.mapping[gen]['field_range']
        if obs is None:
            return [None, None]
        obs = obs.lower()
        m = self.mapping[gen]['mjd_range'][obs]
        if m is None:
            return [None, None]
        return m
        
    def check(self, fieldid = None, mjd = None, obs = None, generations = None):
        """
        Checks if the given fieldid and mjd are valid for the given generations and observatory
        Parameters:
        -----------
        fieldid: int
            the fieldid to check
        mjd: int
            the mjd to check       
        obs: str
            the observatory to check (apo, lco)
        generations: list
            the generations to check (sdssv, as5, legacy, plates, fps, etc)
        Returns:
        --------
        valid: bool
            True if the fieldid and mjd are valid for the given generations and observatory, False otherwise
        """
        if self.mapping is None:
            self.load()
        valid = True
        from boss_drp.Config import config, update_key
        if fieldid is None and mjd is None:
            return False
        if generations is None:
            generations = [g for g in self.mapping.keys() if config.pipe[f'SDSS_Generation.{g}']]
        elif isinstance(generations, str):
            generations = [generations]
        if fieldid is not None:
            if config.pipe['fmjdselect.field'] is not None:
                if isinstance(config.pipe['fmjdselect.field'], (list, tuple)):
                    if int(fieldid) not in [int(f) for f in config.pipe['fmjdselect.field']]:
                        return False
                elif int(config.pipe['fmjdselect.field']) != int(fieldid):
                    return False
            if config.pipe['fmjdselect.fieldrange'] is not None:
                if isinstance(config.pipe['fmjdselect.fieldrange'][0], (list, tuple)):
                    for f in config.pipe['fmjdselect.fieldrange'][0]:
                        if f is None:
                            continue
                        if int(f) > int(fieldid):
                            return False
                elif config.pipe['fmjdselect.fieldrange'][0] is not None:
                    if int(config.pipe['fmjdselect.fieldrange'][0]) > int(fieldid):
                        return False

                if isinstance(config.pipe['fmjdselect.fieldrange'][1], (list, tuple)):
                    for f in config.pipe['fmjdselect.fieldrange'][1]:
                        if f is None:
                            continue
                        if int(f) < int(fieldid):
                            return False
                elif config.pipe['fmjdselect.fieldrange'][1] is not None:
                    if int(config.pipe['fmjdselect.fieldrange'][1]) < int(fieldid):
                        return False

            for gen in generations:
                if config.pipe[f'SDSS_Generation.{gen}']:
                    if self.get(gen, field = True, obs=obs) == [None, None]:
                        continue
                    if fieldid < self.get(gen, field = True, obs=obs)[0]:
                        valid = False
                        continue
                    if fieldid > self.get(gen, field = True, obs=obs)[1]:
                        valid = False
                        continue
                    valid = True
            if not valid:
                return False

        if mjd is not None:
            if config.pipe['fmjdselect.mjd'] is not None:
                if isinstance(config.pipe['fmjdselect.mjd'], list):
                    valid = False
                    for _mjd in config.pipe['fmjdselect.mjd']:
                        if _mjd is None:
                            valid = True
                            continue
                        if int(_mjd) == int(mjd):
                            valid = True
                            break
                    if not valid:
                        return False
                elif int(config.pipe['fmjdselect.mjd']) != int(mjd):
                    return False
            if config.pipe['fmjdselect.mjdrange'] is not None:
                if isinstance(config.pipe['fmjdselect.mjdrange'][0], (list, tuple)):
                    for r in config.pipe['fmjdselect.mjdrange']:
                        if r[0] is None:
                            r[0] = -99999
                        if r[1] is None:
                            r[1] = 1000000000
                        if int(mjd) >= int(r[0]) and int(mjd) <= int(r[1]):
                            valid = True
                    if not valid:
                        return False
                elif int(mjd) < int(config.pipe['fmjdselect.mjdrange'][0]) or int(mjd) > int(config.pipe['fmjdselect.mjdrange'][1]):
                    return False

            genvalid = False
            for gen in generations:
                if config.pipe[f'SDSS_Generation.{gen}']:
                    mjdrange = self.get(gen, mjd = True, obs=obs)
                    if mjdrange == [None, None]:
                        genvalid = True
                        continue
                    if mjdrange[0] is None: 
                        mjdrange[0] = 0
                    if mjdrange[1] is None: 
                        mjdrange[1] = 1000000000
                    if mjdrange[0] <= mjd <= mjdrange[1]:
                        genvalid = True
            if not genvalid:
                valid = False
        return valid

    def find_generation(self, fieldid = None, mjd = None, obs = None):
        """
        Finds the generation for the given fieldid and mjd
        Parameters:
        -----------
        fieldid: int
            the fieldid to find the generation for
        mjd: int
            the mjd to find the generation for       
        obs: str
            the observatory to find the generation for (apo, lco)
        Returns:
        --------
        gen: str
            the generation for the given fieldid and mjd, None if not found
        """
        if self.mapping is None:
            self.load()
        if fieldid is None and mjd is None:
            return None
        gens = []
        for gen in self.mapping.keys():
            if fieldid is not None:
                fieldid = int(fieldid)
                if self.get(gen, field = True, obs=obs) == [None, None]:
                    continue
                if fieldid < self.get(gen, field = True, obs=obs)[0]:
                    continue
                if fieldid > self.get(gen, field = True, obs=obs)[1]:
                    continue
                gens.append(gen)
            if mjd is not None:
                mjd = int(mjd)
                if self.get(gen, mjd = True, obs=obs) == [None, None]:
                    continue
                if mjd < self.get(gen, mjd = True, obs=obs)[0]:
                    continue
                if mjd > self.get(gen, mjd = True, obs=obs)[1]:
                    continue
                gens.append(gen)
        if len(gens) == 0:
            return None
        gens = list(set(gens))
        return gens

    def set_config(self, obs, field_only=False):
        """
        Sets the config.pipe['fmjdselect.field'] and config.pipe['fmjdselect.mjd'] to the appropriate values for the given observatory
        Parameters:
        -----------
        obs: str
            the observatory to set the config for (apo, lco)
        field_only: bool
            if True, only set the field range, not the mjd range
        """
        #TODO: Deal with SDSSV split MJD ends APO and LCO
        from boss_drp.Config import config, update_key
        if self.mapping is None:
            self.load()
        if isinstance(obs, (list, tuple, set)):
            if len(obs) > 1:
                obs = 'apo'
            else:
                obs = obs[0]
        if config.pipe['SDSS_Generation.sdssv']:
            if obs.lower() == 'apo':
                update_key(config.pipe, 'plates', True)
                update_key(config.pipe, 'fps', True)

        gens = []
        if config.pipe['SDSS_Generation.sdssv']:
            gens.append('sdssv')
        if config.pipe['SDSS_Generation.as5']:
            gens.append('as5')
        if config.pipe['SDSS_Generation.legacy']:
            gens.append('legacy')
        if len(gens) == 0:
            gens = ['plates', 'fps']
        for gen in gens:
            if config.pipe[f'SDSS_Generation.{gen}']:
                if (config.pipe['fmjdselect.field'] is None):
                    if config.pipe['fmjdselect.field'] is None:
                        if (config.pipe['fmjdselect.fieldrange'] is None):
                            update_key(config.pipe, 'fieldrange', [self.get(gen, field=True)])
                        else:
                            if not isinstance(config.pipe['fmjdselect.fieldrange'][0], (list, tuple)):
                                update_key(config.pipe, 'fieldrange', [config.pipe['fmjdselect.fieldrange']])
                            if self.get(gen, field=True) != [None, None]:
                                update_key(config.pipe, 'fieldrange', config.pipe['fmjdselect.fieldrange'] + [self.get(gen, field=True)])
                elif not isinstance(config.pipe['fmjdselect.field'], (list, tuple)):
                    config.pipe['fmjdselect.field'] = [config.pipe['fmjdselect.field']]

                if (config.pipe['fmjdselect.mjd'] is None):
                    if config.pipe['fmjdselect.mjd'] is None:
                        if (config.pipe['fmjdselect.mjdrange'] is None):
                            update_key(config.pipe, 'mjdrange', [self.get(gen, mjd=True)])
                        else:
                            if not isinstance(config.pipe['fmjdselect.mjdrange'][0], (list, tuple)):
                                update_key(config.pipe, 'mjdrange', [config.pipe['fmjdselect.mjdrange']])
                            if self.get(gen, mjd=True) != [None, None]:
                                update_key(config.pipe, 'mjdrange', config.pipe['fmjdselect.mjdrange'] + [self.get(gen, mjd=True)])
                    elif not isinstance(config.pipe['fmjdselect.mjd'], (list, tuple)):
                        config.pipe['fmjdselect.mjd'] = [config.pipe['fmjdselect.mjd']]


generations = Generations()