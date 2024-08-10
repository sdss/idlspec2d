import astropy.time
import datetime


class ConvertibleValue:
    def __init__(self, value):
        self.value = value

    def astype(self, totype):
        if not isinstance(totype, type):
            raise ValueError("totype must be a type")
        # Convert to the desired type
        if totype == str:
            return str(int(self.value))  # Convert to int first, then to str
        return totype(self.value)

    def __float__(self):
        return float(self.value)

    def __int__(self):
        return int(self.value)

    def __repr__(self):
        return str(self.value)

class JDATE:
    def __init__(self):
        try:
            self.mjd = astropy.time.Time(datetime.datetime.now(datetime.UTC)).jd - 2400000.5
        except:
            self.mjd = float(astropy.time.Time(datetime.datetime.utcnow()).jd) - 2400000.5
        self._apo = ConvertibleValue(self.sjd(offset=.3))
        self._lco = ConvertibleValue(self.sjd(offset=.4))

    def sjd(self, offset=0.3):
        return self.mjd + offset

    @property
    def apo(self):
        return self._apo

    @property
    def lco(self):
        return self._lco

    def astype(self, totype):
        if not isinstance(totype, type):
            raise ValueError("totype must be a type")
        if totype == str:
            return str(int(self.mjd))
        return totype(self.mjd)

# Example usage
jdate = JDATE()
