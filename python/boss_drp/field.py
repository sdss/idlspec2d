#!/usr/bin/env python3

def field_to_string(field):
    return(str(field).zfill(6))

class Fieldtype:
    def __init__(self, fieldid=None, mjd=None):
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
        
        if fieldid is not None:
            if fieldid == 0:
                self.fps=True
                self.engineering=True
                self.bad=True
            elif int(fieldid) < 15000:
                self.legacy=True
            elif int(fieldid) < 16000:
                self.plates=True
            elif int(fieldid) < 100000:
                self.commissioning = self.fps = True
            else:
                self.fps=True
        elif mjd is not None:
            if int(mjd) == -1:
                self.fps=True
                self.engineering=True
                self.bad=True
            elif int(mjd) < 59030:
                self.legacy=True
            elif int(mjd) < 59550:
                self.plates=True
            else:
                self.fps=True

        if self.fps:
            self.mjd_range=[59550,70000]
            self.field_range=[16000,999999]
        elif self.legacy:
            self.mjd_range=[0,59030]
            self.field_range=[1,14999]
        elif self.plates:
            self.mjd_range=[59030,59550]
            self.field_range=[15000,15999]
        
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
