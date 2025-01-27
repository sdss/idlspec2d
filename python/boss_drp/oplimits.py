import boss_drp

from pydl.pydlutils.yanny import read_table_yanny
import fnmatch
import os
import numpy as np
import copy
def color2hex(color):
    color = color.strip().upper()
    if color == 'RED':
        return '#FF0000'
    if color == 'YELLOW':
        return '#909000'
    return 'black'

def is_number(string):
    try:
        float(string)  # Try converting to a float
        return True
    except ValueError:
        return False

def wildcard_filter(array, pattern):
    """
    Filters an array using fnmatch, handling wildcard matching in both the array
    and the input pattern. Elements in the array like '*' or 'b*' are treated as wildcards.
    """
    if pattern == '*':
        return np.ones(len(array), dtype=bool)  # Match everything if pattern is '*'
    
    # Match pattern against each element, treating '*' and other wildcards in the array
    return np.array([
        fnmatch.fnmatch(pattern, item) or fnmatch.fnmatch(item, pattern) for item in array
    ])

class Limit_query:
    def __init__(self, flavor, field, camera, value, html=False, format = None):
        self.flavor = flavor.lower()
        self.field = field.lower()
        self.camera = camera.lower()
        self.value = value
        self.html = html
        self.format = format

    def __deepcopy__(self, memo):
        # Create a new instance
        new_instance = self.__class__.__new__(self.__class__)
        # Deep copy all attributes
        memo[id(self)] = new_instance  # Add to memo to handle circular references
        for key, value in vars(self).items():
            setattr(new_instance, key, copy.deepcopy(value, memo))
        return new_instance
        
    def copy(self):
        # Use the deepcopy method to return a new copy
        return copy.deepcopy(self)

class opLimits:
    fpath = os.path.join(boss_drp.idlspec2d_dir,'examples','opLimits.par')

    
    def __init__(self):
        self.numlimits = self.norm_table(read_table_yanny(self.fpath,'SPECLIMIT'))
        self.textlimits = self.norm_table(read_table_yanny(self.fpath,'TEXTLIMIT'))
        self.limits = None
        self._limits = None
        self.filter_val = None
        self._arr = False

    def norm_table(self, table):
        for col in table.colnames:
            if table[col].dtype.kind in {'U', 'S'}:  # Check for string columns
                table[col] = [s.lower() for s in table[col]]
        return(table)


    def filter(self, flavor, field, camera):
        field = field.lower().strip()
        flavor = flavor.lower().strip()
        camera = camera.lower().strip()
        self.limits = self.limits[(
                    wildcard_filter(self.limits['field'], field) &
                    wildcard_filter(self.limits['flavor'], flavor) &
                    wildcard_filter(self.limits['camera'], camera)
                    )]
        self._limits = self.limits.copy() if self._arr else None

        
    def filter_val_num(self, value):
        value = float(value)
        self.limits = self.limits[((self.limits['lovalue'] <= value) &  # Value within range
                                   (value <= self.limits['hivalue']))]  # Value within range
    def filter_val_str(self, value):
        value = str(value).lower().strip()
        self.limits = self.limits[wildcard_filter(self.limits['strval'], value)]
      
    def float2str(self, value, format):
        fval = value
        if format is not None:
            if 'd' in format:
                fval = int(fval)
            fval = format.format(fval).strip()

        else:
            fval = str(fval).strip()
        if fval == 'nan':
            fval = '-'
        return(fval)

    def get_type(self, flavor, field, camera):
        self.limits = self.textlimits.copy
        self.filter(flavor,field,camera)
        type = 'num' if len(self.limits) == 0 else 'text'
        self.limits = None
        self._limits = None
        return(type)

    def check(self, flavor=None, field=None, camera=None, value=None, ql=None,
              html=False, type = None, format=None):
        self.limits = None
        if ql is None:
            ql = Limit_query(flavor, field, camera, value, html=html, format = format)
        if ql.value is None or (isinstance(ql.value, (list, tuple)) and len(ql.value) == 0):
            if ql.html:
                return '<span></span>'
            self.limits = None
            return ''
        elif (isinstance(ql.value, (list, tuple)) and len(ql.value) > 0):
            results = []
            if type is None:
                value = str(ql.value[0])
                type = 'text' if not is_number(value) else 'num'
            self._arr = True
            for v in ql.value:
                ql1 = ql.copy()
                ql1.value = v
                print(self._limits, type)
                results.append(self.check(ql = ql1, type = type))
            self._arr = False
            self._limits = None
            return(results)

        if type is None:
            value = str(ql.value)
            type = 'text' if not is_number(value) else 'num'
            
        if type == 'text':
            self.limits = self.textlimits.copy()
            self.filter_val = self.filter_val_str
        else:
            self.limits = self.numlimits.copy()
            self.filter_val = self.filter_val_num

        if self._limits is not None:
            self.limits = self._limits.copy()
        if self.limits is not None:
            self.filter(ql.flavor, ql.field, ql.camera)
            
        self.filter_val(ql.value)
        
        markstr = ''
        if ql.html:
            markstr = '<span> '+self.float2str(ql.value, ql.format)+' </span>'

        if len(self.limits) > 0:
            if len(self.limits) > 1:
                if len(self.limits[self.limits['color'] == 'red']) >= 1 :
                    self.limits = self.limits[self.limits['color'] == 'red'][0]
                else:
                    self.limits = self.limits[0]
            else:
                self.limits = self.limits[0]
            markstr = self.limits['color']
            if ql.html:
                fval = self.float2str(ql.value, ql.format)
                markstr = f'<span style="color:{color2hex(markstr)};font-weight:bold;">'+fval+'</span>'
        self.limits = None
        if not self._arr:
            self._limits = None
        return markstr


oplimits = opLimits()


