from .generations import generations

class Fieldtype:
    def __init__(self, fieldid=None, mjd=None, obs=None):
        self.fieldid=fieldid
        self.mjd=mjd
        self.legacy=False
        self.plates=False
        self.fps=False
        self.dither=False
        self.commissioning=False
        self.bad=False
        self.engineering=False
        self.mjd_range=None
        self.field_range=None
        self.rm_plate = False
        self.obs = obs
        
        if fieldid is not None:
            if fieldid == 0:
                self.fps=True
                self.engineering=True
                self.bad=True
            elif int(fieldid) < generations.get('legacy', field = True, obs=obs)[1]:
                self.legacy=True
            elif int(fieldid) < generations.get('plates', field = True, obs=obs)[1]:
                self.plates=True
            elif int(fieldid) < generations.get('commissioning', field = True, obs=obs)[1]:
                self.commissioning = self.fps = True
            else:
                self.fps=True
        elif mjd is not None:
            if int(mjd) == -1:
                self.fps=True
                self.engineering=True
                self.bad=True
            elif int(mjd) < (generations.get('legacy', mjd=True, obs=self.obs)[1] or 0): 
                self.legacy=True
            elif int(mjd) < (generations.get('plates', mjd=True, obs=self.obs)[1] or 0):
                self.plates=True
            else:
                self.fps=True

        if self.fps:
            self.mjd_range=generations.get('fps', mjd= True, obs=self.obs)
            self.field_range=generations.get('fps', field= True, obs=self.obs)
        elif self.legacy:
            self.mjd_range=generations.get('legacy', mjd= True, obs=self.obs)
            self.field_range=generations.get('legacy', field= True, obs=self.obs)
        elif self.plates:
            self.mjd_range=generations.get('plates', mjd= True, obs=self.obs)
            self.field_range=generations.get('plates', field= True, obs=self.obs)
        
        if self.fieldid in [20903,20931, 20933, 20939, 20955, 20957,20959,20963, 20965,20971, 20973, 20979, 20981,20987, 20989, 21310, 21324, 21325, 22744,22746]:
            self.dither = True

        if self.fieldid in [15000,15001,15002,15038,15070,15071,15171,15172,15173,15252,15253]:
            self.rm_plate = True
        #Define bad
            #no boss fibers
            #if field == 16174: types.bad=True
            #Incorrect design/configuration
            #if field == 16165 & mjd == 59615: types.bad=True
            #unguided
            #if field == 20549 & mjd == 59623: types.bad=True
    def string(self):
        fstr = []

        if self.legacy:
            fstr.append('legacy')
        if self.plates:
            fstr.append('plates')
        if self.fps:
            fstr.append('fps')
        if self.dither:
            fstr.append('dither')
        if self.commissioning:
            fstr.append('commissioning')
        if self.bad:
            fstr.append('bad')
        if self.engineering:
            fstr.append('engineering')
        return ','.join(fstr)
        
    
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        fstr = []

        if self.legacy:
            fstr.append('legacy')
        if self.plates:
            fstr.append('plates')
        if self.fps:
            fstr.append('fps')
        if self.dither:
            fstr.append('dither')
        if self.commissioning:
            fstr.append('commissioning')
        if self.bad:
            fstr.append('bad')
        if self.engineering:
            fstr.append('engineering')
        return ','.join(fstr)
