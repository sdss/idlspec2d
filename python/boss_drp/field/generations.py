import boss_drp
import json
from pathlib import Path

class Generations:
    def __init__(self):
        self.mapping = None
        self._map_file = Path(boss_drp.__file__).parent / 'etc' / 'generations.json'
    def load(self):
        if self.mapping is None:
            with open(self._map_file , "r") as f:
                self.mapping = json.load(f)

    def get(self, gen, field = False, mjd = False, obs= None):
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
        
    def set_config(self, obs, field_only=False):
        from boss_drp.Config import config, update_key
        if self.mapping is None:
            self.load()
        if isinstance(obs, (list, tuple, set)):
            if len(obs) > 1:
                obs = 'apo'
            else:
                obs = obs[0]
        for gen in self.mapping.keys():
            if config.pipe[f'SDSS_Generation.{gen}']:
                if (config.pipe['fmjdselect.field'] is None):
                    if (config.pipe['fmjdselect.fieldstart'] is None):
                        update_key(config.pipe, 'fieldstart', self.mapping[gen]['field_range'][0])
                    if (config.pipe['fmjdselect.fieldend'] is None):
                        update_key(config.pipe, 'fieldend', self.mapping[gen]['field_range'][1])
            if field_only:
                continue
            if (config.pipe['fmjdselect.mjd'] is None):
                if obs is None:
                    return
                obs = obs.lower()
                if self.mapping[gen]['mjd_range'][obs] is None:
                    return
                if (config.pipe['fmjdselect.mjdstart'] is None):
                    update_key(config.pipe, 'mjdstart', self.mapping[gen]['mjd_range'][obs][0])
                if (config.pipe['fmjdselect.mjdend'] is None):
                    update_key(config.pipe, 'mjdend', self.mapping[gen]['mjd_range'][obs][1])

generations = Generations()